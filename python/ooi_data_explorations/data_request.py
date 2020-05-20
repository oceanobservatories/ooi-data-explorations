#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
import pytz
import xarray as xr
import dateutil.parser as parser

from ooi_data_explorations import M2M_URLS
from ooi_data_explorations.common import dr_inputs, m2m_request, m2m_collect, get_deployment_dates


def filter_urls(site, assembly, instrument, method):
    """
    Takes the M2M_URLS dictionary and searches for the instrument of interest
    as defined by the site code, assembly type, instrument class, and data
    delivery method to return the OOI specific site, node and stream names
    needed to request the data.

    :param site: OOI eight letter site code (e.g. CE04OSPS for the Oregon
        Offshore Shallow Profiler)
    :param assembly: Assembly grouping name (e.g. midwater for the 200 m
        Platform)
    :param instrument: The instrument class name (e.g. phsen for the
        SAMI2-pH sensor)
    :param method: The data delivery method (e.g. streamed for cabled
        streaming data)

    :return node: The OOI specific node code(s) for the assembly
    :return sensor: The OOI specific sensor code(s) for the instrument class
    :return stream: The OOI specific stream name(s) for the site, node, sensor
        and delivery method combination
    """
    # Use the site, assembly, instrument, and data delivery method to determine the OOI specific site, node, sensor,
    # method and stream
    node = []
    sensor = []
    stream = []

    # pare the larger dictionary down to the site of interest and check if a valid site was used
    m2m_urls = M2M_URLS.get(site.upper())
    if not m2m_urls:
        raise SyntaxError('Unknown site code: %s' % site)

    # make sure the correct data delivery method was specified
    if method not in ['streamed', 'telemetered', 'recovered_host', 'recovered_inst', 'recovered_cspp', 'recovered_wfp']:
        raise SyntaxError('Unknown data delivery method: %s' % method)

    # find the instrument(s) of interest in the assembly group
    for grouping in m2m_urls.get('assembly'):
        if grouping.get('type') == assembly or grouping.get('subassembly') == assembly:
            for instrmt in grouping['instrument']:
                if instrmt['class'] == instrument:
                    node.append(instrmt.node)
                    sensor.append(instrmt.sensor)
                    stream.append(instrmt.stream.get(method))

    # check to see if we were able to find the system of interest
    if not stream:
        raise RuntimeWarning('Instrument defined by %s-%s-%s-%s cannot be found.' % (site, assembly, instrument,
                                                                                     method))

    # return the OOI specific names for the node(s), sensor(s) and stream(s)
    return node, sensor, stream


def data_request(site, assembly, instrument, method, **kwargs):
    """
    Requests data via the OOI M2M API using the site code, assembly type,
    instrument class and data delivery method as defined in the m2m_urls.yml to
    construct the OOI specific data request.

    :param site: OOI site code as an 8 character string
    :param assembly: The assembly type where the instrument is located
    :param instrument: the OOI instrument class name for the instrument of
        interest
    :param method: The data delivery method for the system of interest
    :param kwargs: Takes the following optional keyword arguments:

        start: Starting date/time for the data request in a dateutil.parser
            recognizable form. If ``None``, the default, the beginning of the
            data record will be used
        stop: Ending date/time for the data request in a dateutil.parser
            recognizable form If ``None``, the default, the end of the data
            record will be used
        deploy: Use the deployment number, entered as an integer, to set the
            starting and ending dates If ``None``, the default, the starting
            and ending dates are used. If you enter both, the deployment number
            will take priority in setting the start and end dates
        aggregate: In cases where more than one instance of an instrument class
            is part of an assembly, will collect all of the data if the integer
            value entered is ``0``, or the specific instance of the instrument
            is requested if any value greater than ``0`` is used. If ``None``,
            the default, the first instance of an instrument will be used.

    :return data: Returns the request data as an xarray dataset for further analysis
    """
    # setup inputs to the function, make sure case is correct
    site = site.upper()
    assembly = assembly.lower()
    instrument = instrument.lower()
    method = method.lower()

    # parse the keyword arguments
    start = None
    stop = None
    deploy = None
    aggregate = None
    for key, value in kwargs.items():
        if key not in ['start', 'stop', 'deploy', 'aggregate']:
            raise KeyError('Unknown keyword (%s) argument.' % key)
        else:
            if key == 'start':
                start = value
            if key == 'stop':
                stop = value
            if key == 'deploy':
                deploy = value
            if key == 'aggregate':
                aggregate = value

    # use the assembly, instrument and data delivery methods to find the system of interest
    node, sensor, stream = filter_urls(site, assembly, instrument, method)

    # check the formatting of the start and end dates. We need to be able to parse and convert to an ISO format.
    if start:
        # check the formatting of the start date string and convert to the ISO format used by the M2M API
        try:
            start = parser.parse(start)
            start = start.astimezone(pytz.utc)
            start = start.strftime('%Y-%m-%dT%H:%M:%S.000Z')
        except parser.ParserError:
            raise SyntaxError('Formatting of the starting date string needs to be in a recognizable format')

    if stop:
        # check the formatting of the stop date string and convert to the ISO format used by the M2M API
        try:
            stop = parser.parse(stop)
            stop = stop.astimezone(pytz.utc)
            stop = stop.strftime('%Y-%m-%dT%H:%M:%S.000Z')
        except parser.ParserError:
            raise SyntaxError('Formatting of the ending date string needs to be in a recognizable format')

    if deploy:
        # Determine start and end dates based on the deployment number
        start, stop = get_deployment_dates(site, node[0], sensor[0], deploy)
        if not start or not stop:
            exit_text = ('Deployment dates are unavailable for %s-%s-%s-%s, deployment %02d.' % (site.lower(), assembly,
                                                                                                 instrument, method,
                                                                                                 deploy))
            raise RuntimeWarning(exit_text)

    # for some cases, there maybe more than 1 stream, but in general, we only want the first one
    if isinstance(stream[0], list):
        stream = stream[0][0]
    else:
        stream = stream[0]

    tag = ('.*{instrument}.*\\.nc$'.format(instrument=instrument.upper()))  # set regex tag to use when downloading
    data = None     # setup the default data set

    # check if there are multiple instances of this instrument class on the assembly
    if len(node) > 1:
        print('There are multiple instances of the instrument %s under %s-%s.' % (instrument, site.lower(), assembly))

    # check if we are aggregating the multiple instruments into a single data set
    if isinstance(aggregate, int):
        if aggregate == 0:
            # request all of the instruments associated with this site, assembly, instrument and method
            print(('Requesting all %d instances of this instrument. Data sets will be concatenated\n'
                   'and a new variable called `sensor_count` will be added to help distinguish the \n'
                   'instruments for later processing.') % len(node))
            for i in range(len(node)):
                r = m2m_request(site, node[i], sensor[i], method, stream, start, stop)
                if r:
                    temp = m2m_collect(r, tag)
                    temp['sensor_count'] = temp['deployment'] * 0 + i + 1
                    if not data:
                        data = temp
                    else:
                        data = xr.concat([data, temp], dim='time')
        else:
            # request a specific instrument of the multiple instruments associated with this site, assembly,
            # instrument and method.
            if aggregate > len(node):
                raise SyntaxError('Only %d instruments available, you selected %d' % (len(node), aggregate))

            print('Requesting instrument %d out of %d.' % (aggregate, len(node)))
            i = aggregate - 1
            r = m2m_request(site, node[i], sensor[i], method, stream, start, stop)
            if r:
                data = m2m_collect(r, tag)
    else:
        r = m2m_request(site, node[0], sensor[0], method, stream, start, stop)
        if r:
            data = m2m_collect(r, tag)

    if not data:
        raise RuntimeWarning('Data unavailable for %s-%s-%s-%s.' % (site.lower(), assembly, instrument, method))

    # return the resulting data, which is an xarray.Dataset object
    for v in data.variables:
        # first convert strings with data types set as objects or S64 with binary encoding
        if data[v].dtype == np.dtype('O') or data[v].dtype == np.dtype('S64'):
            data[v] = data[v].astype(np.str)

    return data


def main(argv=None):
    """
    CLI interface to the data_request function. Requests data via the OOI M2M
    API using the site code, assembly type, instrument class and data delivery
    method as defined in the m2m_urls.yml to construct the OOI specific data
    request.

    :param argv: Command line inputs used by the function are:

        site (-s, --site): OOI site code as an 8 character string (required)
        assembly (-a, --assembly): Assembly or subassembly type where the
            instrument is located (required)
        instrument (-i, --instrument): OOI instrument class name for the
            instrument of interest (required)
        method (-m, --method): The data delivery method for the system of
            interest (required)
        outfile (-o, --outfile): An absolute or relative path and file name
            for where the resulting data should be saved as a NetCDF file
            (required)

        start (-bt, --beginDT): Starting or beginning date/time for the data
            request in a dateutil.parser recognizable form. If ``None``, the
            default, the beginning of the data record will be used (optional)
        stop (-et, --endDT): Ending date/time for the data request in a
            dateutil.parser recognizable form If ``None``, the default, the end
            date of the data record will be used (optional)
        deploy (-dp, --deploy): Use the deployment number, entered as an
            integer, to set the beginning and ending dates If ``None``, the
            default, the beginning and ending dates are used. If you enter
            both, the deployment number will take priority in setting the start
            and end dates (optional)
        aggregate (-ag, --aggregate): In cases where more than one instance of
            an instrument class is part of an assembly, collect all of the data
            if the integer value entered is ``0``, or the specific instance of
            the instrument is requested if any value greater than ``0`` is
            used. If ``None``, the default, the first instance of an instrument
            will be requested. (optional)

    :return: Saves the data to disk as a NetCDF file.
    """
    args = dr_inputs(argv)
    site = args.site
    assembly = args.assembly
    instrument = args.instrument
    method = args.method
    outfile = os.path.abspath(args.outfile)
    deploy = args.deploy
    start = args.start
    stop = args.stop
    aggregate = args.aggregate

    # request the data
    data = data_request(site, assembly, instrument, method, start=start, stop=stop, deploy=deploy, aggregate=aggregate)

    # save the data to disk
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    data.to_netcdf(outfile, mode='w', format='NETCDF4', engine='h5netcdf')


if __name__ == '__main__':
    main()
