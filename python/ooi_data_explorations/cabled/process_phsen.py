#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import xarray as xr

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, get_deployment_dates, get_vocabulary, \
    dt64_epoch, update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_phsen import PHSEN


def phsen_streamed(ds):
    """
    Takes PHSEN data streamed from instruments deployed by the Regional Cabled Array and cleans up the data set to make
    it more user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that
    are of no use/value. Additionally, re-organize some of the variables to permit better assessments of the data.

    :param ds: initial PHSEN data set recorded by the data logger system and downloaded from OOI via the M2M system
    :return: cleaned up and reorganized data set
    """
    # drop some of the variables:
    #   checksum == meaningless
    #   record_type == there is only one, don't need this
    #   record_length == meaningless
    #   signal_intensity_434, part of the light measurements array, redundant so can remove
    #   signal_intensity_578, part of the light measurements array, redundant so can remove
    #   provenance == better to access with direct call to OOI M2M api, it doesn't work well in this format
    ds = ds.reset_coords()
    ds = ds.drop(['checksum', 'record_type', 'record_length', 'signal_intensity_434',
                  'signal_intensity_578', 'provenance'])

    # convert the internal_timestamp values from a datetime64[ns] object to a floating point number with the time in
    # seconds, replacing the internal_timestamp with the record_time (the internal_timestamp is incorrectly set in the
    # NetCDF file).
    ds['internal_timestamp'] = ('time', dt64_epoch(ds.record_time))
    ds['internal_timestamp'].attrs = dict({
        'long_name': 'Internal SAMI-pH Clock Time',
        'standard_name': 'time',
        'units': 'seconds since 1970-01-01 00:00:00 0:00',
        'calendar': 'gregorian',
        'comment': ('Comparing the instrument internal clock versus the GPS referenced sampling time will allow for ' +
                    'calculations of the instrument clock offset and drift. Useful when working with the ' +
                    'recovered instrument data where no external GPS referenced clock is available.')
    })
    ds = ds.drop(['record_time'])

    # rename some of the variables for better clarity
    rename = {
        'voltage_battery': 'raw_battery_voltage',
        'thermistor_start': 'raw_thermistor_start',
        'thermistor_end': 'raw_thermistor_end',
        'phsen_thermistor_temperature': 'thermistor_temperature',
        'phsen_battery_volts': 'battery_voltage',
        'ph_seawater': 'seawater_ph',
        'ph_seawater_qc_executed': 'seawater_ph_qc_executed',
        'ph_seawater_qc_results': 'seawater_ph_qc_results'
    }
    ds = ds.rename(rename)

    # now we need to reset the light and reference arrays to named variables that will be more meaningful and useful in
    # the final data files
    nrec = len(ds['time'].values)
    light = np.array(np.vstack(ds['ph_light_measurements'].values), dtype='int32')
    light = np.atleast_3d(light)
    light = np.reshape(light, (nrec, 23, 4))  # 4 sets of 23 seawater measurements
    reference_434 = light[:, :, 0]            # reference signal, 434 nm
    signal_434 = light[:, :, 1]               # signal intensity, 434 nm (PH434SI_L0)
    reference_578 = light[:, :, 2]            # reference signal, 578 nm
    signal_578 = light[:, :, 3]               # signal intensity, 578 nm (PH578SI_L0)

    refnc = np.array(np.vstack(ds['reference_light_measurements'].values), dtype='int32')
    refnc = np.atleast_3d(refnc)
    refnc = np.reshape(refnc, (nrec, 4, 4))   # 4 sets of 4 DI water measurements (blanks)
    blank_refrnc_434 = refnc[:, :, 0]  # DI blank reference, 434 nm
    blank_signal_434 = refnc[:, :, 1]  # DI blank signal, 434 nm
    blank_refrnc_578 = refnc[:, :, 2]  # DI blank reference, 578 nm
    blank_signal_578 = refnc[:, :, 3]  # DI blank signal, 578 nm

    # create a data set with the reference and light measurements
    ph = xr.Dataset({
        'blank_refrnc_434': (['time', 'blanks'], blank_refrnc_434.astype('int32')),
        'blank_signal_434': (['time', 'blanks'], blank_signal_434.astype('int32')),
        'blank_refrnc_578': (['time', 'blanks'], blank_refrnc_578.astype('int32')),
        'blank_signal_578': (['time', 'blanks'], blank_signal_578.astype('int32')),
        'reference_434': (['time', 'measurements'], reference_434.astype('int32')),
        'signal_434': (['time', 'measurements'], signal_434.astype('int32')),
        'reference_578': (['time', 'measurements'], reference_578.astype('int32')),
        'signal_578': (['time', 'measurements'], signal_578.astype('int32'))
    }, coords={'time': ds['time'], 'measurements': np.arange(0, 23).astype('int32'),
               'blanks': np.arange(0, 4).astype('int32')
               })
    ds = ds.drop(['light_measurements', 'reference_light_measurements'])

    # merge the data sets back together
    ds = ds.merge(ph)

    # reset some of the variable attributes, and ...
    for v in ds.variables:  # variable level attributes
        if v in PHSEN:
            ds[v].attrs = PHSEN[v]

    # ... add the renamed information
    for key, value in rename.items():
        ds[value].attrs['ooinet_variable_name'] = key

    # and reset some of the data types
    data_types = ['deployment', 'raw_thermistor_end', 'raw_thermistor_start', 'unique_id', 'raw_battery_voltage']
    for v in data_types:
        ds[v] = ds[v].astype('int32')

    return ds


def main(argv=None):
    # setup the input arguments
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    method = args.method
    stream = args.stream
    deploy = args.deploy
    start = args.start
    stop = args.stop

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

    # Request the data
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    if not r:
        exit_text = ('Data unavailable for %s-%s-%s, deployment %02d. Check request.' % (site, node, sensor, deploy))
        raise SystemExit(exit_text)

    # Valid request, start downloading the data
    if deploy:
        phsen = m2m_collect(r, ('.*deployment%04d.*PHSEN.*\\.nc$' % deploy))
    else:
        phsen = m2m_collect(r, '.*PHSEN.*\\.nc$')

    if not phsen:
        exit_text = ('Data unavailable for %s-%s-%s. Check request.' % (site, node, sensor))
        raise SystemExit(exit_text)

    # clean-up and reorganize
    phsen = phsen_streamed(phsen)

    vocab = get_vocabulary(site, node, sensor)[0]
    phsen = update_dataset(phsen, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(os.path.join(CONFIG['base_dir']['m2m_base'], args.outfile))
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    phsen.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
