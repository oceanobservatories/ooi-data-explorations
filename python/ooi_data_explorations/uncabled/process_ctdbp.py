#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, get_deployment_dates, get_vocabulary, dt64_epoch, \
    update_dataset, CONFIG, ENCODINGS


def ctdbp_datalogger(ds, burst=False):
    """
    Takes ctdbp data recorded by the data loggers used in the CGSN/EA moorings and cleans up the data set to make
    it more user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that
    are of no use/value.

    :param ds: initial ctdbp data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return: cleaned up data set
    """
    # drop some of the variables:
    #   dcl_controller_timestamp == time, redundant so can remove
    #   date_time_string = internal_timestamp, redundant so can remove
    #   provenance == better to access with direct call to OOI M2M api, it doesn't work well in this format
    ds = ds.reset_coords()
    ds = ds.drop(['dcl_controller_timestamp', 'date_time_string', 'provenance'])

    # convert the time values from a datetime64[ns] object to a floating point number with the time in seconds
    ds['internal_timestamp'] = ('time', dt64_epoch(ds.internal_timestamp))
    ds['internal_timestamp'].attrs = dict({
        'long_name': 'Internal CTD Clock Time',
        'standard_name': 'time',
        'units': 'seconds since 1970-01-01 00:00:00 0:00',
        'calendar': 'gregorian',
        'comment': ('Comparing the instrument internal clock versus the GPS referenced sampling time will allow for ' +
                    'calculations of the instrument clock offset and drift. Useful when working with the ' +
                    'recovered instrument data where no external GPS referenced clock is available.')
    })

    # rename some of the variables for better clarity
    rename = {
        'temp': 'temperature',
        'temp_qc_executed': 'temperature_qc_executed',
        'temp_qc_results': 'temperature_qc_results'
    }
    ds = ds.rename(rename)
    for key, value in rename.items():
        ds[value].attrs['ooinet_variable_name'] = key

    # correct incorrect units
    ds['temperature'].attrs['units'] = 'degree_Celsius'

    # ancillary_variables attribute set incorrectly (should be a space separated list) for certain variables
    ds['temperature'].attrs['ancillary_variables'] = 'temperature_qc_executed temperature_qc_results'
    ds['conductivity'].attrs['ancillary_variables'] = 'conductivity_qc_executed conductivity_qc_results'
    ds['pressure'].attrs['ancillary_variables'] = 'pressure_qc_executed pressure_qc_results'
    ds['practical_salinity'].attrs['ancillary_variables'] = ('conductivity temperature pressure ' +
                                                             'practical_salinity_qc_executed ' +
                                                             'practical_salinity_qc_results')
    ds['density'].attrs['ancillary_variables'] = ('conductivity temperature pressure lat lon ' +
                                                  'density_qc_executed density_qc_results')

    if burst:   # re-sample the data to a defined time interval using a median average
        # create the burst averaging
        burst = ds
        burst['time'] = burst['time'] - np.timedelta64(450, 's')    # center time windows for 15 minute bursts
        burst = burst.resample(time='15Min', keep_attrs=True, skipna=True).median()
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # reset the attributes...which keep_attrs should do...
        burst.attrs = ds.attrs
        for v in burst.variables:
            burst[v].attrs = ds[v].attrs

        # save the newly average data
        ds = burst

    return ds


def ctdbp_instrument(ds, burst=False):
    """
    Takes ctdbp data recorded by internally by the instrument, and cleans up the data set to make it more
    user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that are
    of no use/value.

    :param ds: initial ctdbp data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return: cleaned up data set
    """
    """
    Takes ctdbp data recorded by internally by the instrument, and cleans up the data set to make it more
    user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that are
    of no use/value.

    :param ds: initial ctdbp data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return: cleaned up data set
    """
    # drop some of the variables:
    #   ctd_time == time, redundant so can remove
    #   internal_timestamp == time, redundant so can remove
    #   conductivity_qc_*, pressure_qc_* == raw measurements, no QC tests should be run
    ds = ds.reset_coords()
    ds = ds.drop(['ctd_time', 'internal_timestamp', 'provenance', 'conductivity_qc_executed', 'conductivity_qc_results',
                  'pressure_qc_executed', 'pressure_qc_results'])

    # rename some of the variables for better clarity, two blocks to keep from stepping on ourselves
    rename1 = {
        'temperature': 'raw_temperature',
        'conductivity': 'raw_conductivity',
        'pressure': 'raw_pressure',
        'pressure_temp': 'raw_pressure_temperature',
    }
    rename2 = {
        'ctdbp_seawater_conductivity': 'conductivity',
        'ctdbp_seawater_conductivity_qc_executed': 'conductivity_qc_executed',
        'ctdbp_seawater_conductivity_qc_results': 'conductivity_qc_results',
        'ctdbp_seawater_temperature': 'temperature',
        'ctdbp_seawater_temperature_qc_executed': 'temperature_qc_executed',
        'ctdbp_seawater_temperature_qc_results': 'temperature_qc_results',
        'ctdbp_seawater_pressure': 'pressure',
        'ctdbp_seawater_pressure_qc_executed': 'pressure_qc_executed',
        'ctdbp_seawater_pressure_qc_results': 'pressure_qc_results',
    }
    ds = ds.rename(rename1)
    for key, value in rename1.items():
        ds[value].attrs['ooinet_variable_name'] = key
    ds = ds.rename(rename2)
    for key, value in rename2.items():
        ds[value].attrs['ooinet_variable_name'] = key

    # correct incorrect units and attributes
    ds['temperature'].attrs['units'] = 'degree_Celsius'
    ds['raw_conductivity'].attrs['comment'] = ('Raw seawater conductivity measurement recorded internally by the ' +
                                               'instrument in counts')
    ds['raw_conductivity'].attrs['long_name'] = 'Raw Seawater Conductivity'
    ds['raw_temperature'].attrs['comment'] = ('Raw seawater temperature measurement recorded internally by the ' +
                                              'instrument in counts')
    ds['raw_temperature'].attrs['long_name'] = 'Raw Seawater Temperature'
    ds['raw_pressure'].attrs['comment'] = ('Raw seawater pressure measurement recorded internally by the ' +
                                           'instrument in counts')
    ds['raw_pressure'].attrs['long_name'] = 'Raw Seawater Pressure'
    ds['raw_pressure_temperature'].attrs['comment'] = ('Raw pressure sensor thermistor temperature, internal to the ' +
                                                       'sensor and recorded internally by the instrument in counts. ' +
                                                       'Used to convert the raw pressure measurement, compensating ' +
                                                       'for the sensor temperature, to pressure reported in dbar.')
    ds['raw_pressure_temperature'].attrs['long_name'] = 'Raw Seawater Pressure Sensor Temperature'

    # ancillary_variables attribute set incorrectly (should be a space separated list) for certain variables
    ds['temperature'].attrs['ancillary_variables'] = 'raw_temperature temperature_qc_executed temperature_qc_results'
    ds['conductivity'].attrs['ancillary_variables'] = ('raw_conductivity conductivity_qc_executed ' +
                                                       'conductivity_qc_results')
    ds['pressure'].attrs['ancillary_variables'] = ('raw_pressure raw_pressure_temperature ' +
                                                   'pressure_qc_executed pressure_qc_results')
    ds['practical_salinity'].attrs['ancillary_variables'] = ('conductivity temperature pressure ' +
                                                             'practical_salinity_qc_executed ' +
                                                             'practical_salinity_qc_results')
    ds['density'].attrs['ancillary_variables'] = ('conductivity temperature pressure lat lon ' +
                                                  'density_qc_executed density_qc_results')

    if burst:   # re-sample the data to a defined time interval using a median average
        # create the burst averaging
        burst = ds
        burst['time'] = burst['time'] - np.timedelta64(450, 's')    # center time windows for 15 minute bursts
        burst = burst.resample(time='15Min', keep_attrs=True, skipna=True).median()
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # reset the attributes...which keep_attrs should do...
        burst.attrs = ds.attrs
        for v in burst.variables:
            burst[v].attrs = ds[v].attrs

        # save the newly average data
        ds = burst

    return ds


def main(argv=None):
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    method = args.method
    stream = args.stream
    deploy = args.deploy
    start = args.start
    stop = args.stop
    burst = args.burst

    # determine the start and stop times for the data request based on either the deployment number or user entered
    # beginning and ending dates.
    if not deploy or (start and stop):
        return SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')
    else:
        if deploy:
            # Determine start and end dates based on the deployment number
            start, stop = get_deployment_dates(site, node, sensor, deploy)
            if not start or not stop:
                exit_text = ('Deployment dates are unavailable for %s-%s-%s, deployment %02d.' % (site, node, sensor,
                                                                                                  deploy))
                raise SystemExit(exit_text)

    # Request the data for download
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    if not r:
        exit_text = ('Request failed for %s-%s-%s. Check request.' % (site, node, sensor))
        raise SystemExit(exit_text)

    # Valid request, start downloading the data
    if deploy:
        ctdbp = m2m_collect(r, '.*deployment%04d.*CTDBP.*\\.nc$')
    else:
        ctdbp = m2m_collect(r, '.*CTDBP.*\\.nc$')

    if not ctdbp:
        exit_text = ('Data unavailable for %s-%s-%s. Check request.' % (site, node, sensor))
        raise SystemExit(exit_text)

    # clean-up and reorganize
    if method in ['telemetered', 'recovered_host']:
        ctdbp = ctdbp_datalogger(ctdbp, burst)
    else:
        ctdbp = ctdbp_instrument(ctdbp, burst)

    vocab = get_vocabulary(site, node, sensor)[0]
    ctdbp = update_dataset(ctdbp, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(os.path.join(CONFIG['base_dir']['m2m_base'], args.outfile))
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    ctdbp.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
