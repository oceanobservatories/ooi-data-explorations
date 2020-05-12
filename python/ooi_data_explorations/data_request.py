#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
import xarray as xr

from ooi_data_explorations import M2M_URLS
from ooi_data_explorations.common import inputs, m2m_request, m2m_collect, get_deployment_dates


def filter_urls(site, assembly, instrument, method):
    """
    Takes the master M2M dictionary, used to construct URLs for data requests, and searches for the instrument
    of interest as defined by the site code, assembly type, instrument class, and data delivery method to return
    the OOI specific site, node and stream names needed to request the data.

    :param site: OOI eight letter site code (e.g. CE04OSPS for the Oregon Offshore Shallow Profiler)
    :param assembly: Assembly grouping name (e.g. midwater for the 200 m Platform)
    :param instrument: The instrument class name (e.g. phsen for the SAMI2-pH sensor)
    :param method: The data delivery method (e.g. streamed for cabled streaming data)
    :return node: The OOI specific node code(s) for the assembly
    :return sensor: The OOI specific sensor code(s) for the instrument class
    :return stream: The OOI specific stream name(s) for the site, node, sensor and delivery method combination
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

    # make sure the correct assembly grouping was used
    if assembly not in ['surface_buoy', 'midwater', 'nsif', 'riser', 'seafloor', 'auv', 'glider', 'profiler']:
        raise SyntaxError('Unknown assembly type: %s' % assembly)

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

    :param site:
    :param assembly:
    :param instrument:
    :param method:
    :param kwargs:
    :return:
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

    # determine the start and stop times for the data request based on either the deployment number or user entered
    # beginning and ending dates.
    if not deploy or (start and stop):
        raise SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')

    # use the assembly, instrument and data delivery methods to find the system of interest
    node, sensor, stream = filter_urls(site, assembly, instrument, method)

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

    # Request the data for download
    tag = ('.*{instrument}.*\\.nc$'.format(instrument=instrument.upper()))
    data = None
    if aggregate and len(node) > 1:
        print('There are multiple instances of the instrument %s under %s-%s.' % (instrument, site, assembly))
        if aggregate == 0:
            # request all of the instruments associated with this site, assembly, instrument and method
            print(('Requesting all %d instances of this instrument. Data sets will\n' +
                   'be concatenated and a new variable called ''sensor_count''\n' +
                   'will be added to help distinguish the instruments for later\n' +
                   'processing.') % len(node))
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
    """'''

    :param argv:
    :return:
    """
    args = inputs(argv)
    site = args.site
    assembly = args.assembly
    instrument = args.instrument
    method = args.method
    deploy = args.deploy
    start = args.start
    stop = args.stop
    aggregate = args.aggregate
    outfile = args.outfile

    if not deploy or (start and stop):
        return SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')

    # request the data
    if deploy:
        data = data_request(site, assembly, instrument, method, deploy=deploy, aggregate=aggregate)
    else:
        data = data_request(site, assembly, instrument, method, start=start, stop=stop, aggregate=aggregate)

    # save the data to disk
    outfile = os.path.abspath(outfile)
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    data.to_netcdf(outfile, mode='w', format='NETCDF4', engine='h5netcdf')


if __name__ == '__main__':
    main()
