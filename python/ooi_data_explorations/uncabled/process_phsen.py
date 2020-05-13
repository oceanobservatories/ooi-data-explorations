#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import re
import xarray as xr

from instruments import inputs, m2m_collect, m2m_request, get_deployment_dates, get_vocabulary, \
    dt64_epoch, update_dataset, CONFIG, ENCODINGS

# Setup some attributes, used to replace those incorrectly set, or needed after the processing below
PHSEN = {
    'unique_id': {
        'long_name': 'Instrument Unique ID',
        'comment': ('One byte checksum summary of the instrument serial number, name, calibration date and firmware ' +
                    'version serving as a proxy for a unique ID. While identified as the instrument unique ID, it is ' +
                    'possible for more than one instrument to have the same checksum summary. Thus, the uniqueness ' +
                    'of this value should be considered with a certain degree of caution.'),
        # 'units': ''    # deliberately left blank, no units for this value
    },
    'measurements': {
        'long_name': 'Measurements Array',
        'comment': 'Dimensional indexing array created for the light measurements.',
        # 'units': ''    # deliberately left blank, no units for this value
    },
    'blanks': {
        'long_name': 'Blanks Array',
        'comment': 'Dimensional indexing array created for the reference (blanks) measurements.',
        # 'units': ''    # deliberately left blank, no units for this value
    },
    'blank_refrnc_434': {
        'long_name': 'DI Blank Reference Intensity at 434 nm',
        'comment': ('Optical absorbance reference intensity at 434 nm. Measured with deionized water. Reference and ' +
                    'signal intensities range between 0 and 4096. Values should be greater than ~1500. Lower ' +
                    'intensities will result in higher noise in the absorbance and pH measurements. Obtained from ' +
                    'the reference_light_measurements variable in the Data Portal sourced data file.'),
        'units': 'counts',
        '_FillValue': -9999999,
    },
    'blank_signal_434': {
        'long_name': 'DI Blank Signal Intensity at 434 nm',
        'comment': ('Optical absorbance signal intensity at 434 nm. Measured with deionized water. Reference and ' +
                    'signal intensities range between 0 and 4096. Values should be greater than ~1500. Lower ' +
                    'intensities will result in higher noise in the absorbance and pH measurements. Obtained from ' +
                    'the reference_light_measurements variable in the Data Portal sourced data file.'),
        'units': 'counts',
        '_FillValue': -9999999,
    },
    'blank_refrnc_578': {
        'long_name': 'DI Blank Reference Intensity at 578 nm',
        'comment': ('Optical absorbance reference intensity at 578 nm. Measured with deionized water. Reference and ' +
                    'signal intensities range between 0 and 4096. Values should be greater than ~1500. Lower ' +
                    'intensities will result in higher noise in the absorbance and pH measurements. Obtained from ' +
                    'the reference_light_measurements variable in the Data Portal sourced data file.'),
        'units': 'counts',
        '_FillValue': -9999999,
    },
    'blank_signal_578': {
        'long_name': 'DI Blank Signal Intensity at 578 nm',
        'comment': ('Optical absorbance signal intensity at 578 nm. Measured with deionized water. Reference and ' +
                    'signal intensities range between 0 and 4096. Values should be greater than ~1500. Lower ' +
                    'intensities will result in higher noise in the absorbance and pH measurements. Obtained from ' +
                    'the reference_light_measurements variable in the Data Portal sourced data file.'),
        'units': 'counts',
        '_FillValue': -9999999,
    },
    'reference_434': {
        'long_name': 'Reference Intensity at 434 nm',
        'comment': ('Optical absorbance reference intensity at 434 nm. Reference and signal intensities range ' +
                    'between 0 and 4096. Values should be greater than ~1500. Lower intensities will result in ' +
                    'higher noise in the absorbance and pH measurements. Obtained from the light_measurements' +
                    'variable in the Data Portal sourced data file.'),
        'units': 'counts'
    },
    'signal_434': {
        'long_name': 'Signal Intensity at 434 nm',
        'comment': ('Optical absorbance signal intensity at 434 nm. Reference and signal intensities range between 0 ' +
                    'and 4096. Values should be greater than ~1500. Lower intensities will result in higher noise in ' +
                    'the absorbance and pH measurements. Obtained from the light_measurements variable in the Data ' +
                    'Portal sourced data file.'),
        'data_product_identifier': 'PH434SI_L0',
        'units': 'counts'
    },
    'reference_578': {
        'long_name': 'Reference Intensity at 578 nm',
        'comment': ('Optical absorbance reference intensity at 578 nm. Reference and signal intensities range ' +
                    'between 0 and 4096. Values should be greater than ~1500. Lower intensities will result in ' +
                    'higher noise in the absorbance and pH measurements. Obtained from the light_measurements ' +
                    'variable in the Data Portal sourced data file.'),
        'units': 'counts'
    },
    'signal_578': {
        'long_name': 'Signal Intensity at 578 nm',
        'comment': ('Optical absorbance signal intensity at 578 nm. Reference and signal intensities range between 0 ' +
                    'and 4096. Values should be greater than ~1500. Lower intensities will result in higher noise in ' +
                    'the absorbance and pH measurements. Obtained from the light_measurements variable in the Data ' +
                    'Portal sourced data file.'),
        'data_product_identifier': 'PH578SI_L0',
        'units': 'counts'
    },
    'raw_thermistor_start': {
        'long_name': 'Raw Thermistor Temperature, Measurement Start',
        'comment': 'Thermistor resistivity measured in counts at the start of the measurement cycle.',
        'units': 'counts'
    },
    'raw_thermistor_end': {
        'long_name': 'Raw Thermistor Temperature, Measurement End',
        'comment': 'Thermistor resistivity measured in counts at the end of the measurement cycle.',
        'units': 'counts'
    },
    'raw_battery_voltage': {
        'long_name': 'Raw Battery Voltage',
        'comment': ('Raw internal battery voltage measured in counts. May actually reflect external voltage if ' +
                    'external power is applied'),
        'units': 'counts'
    },
    'thermistor_temperature': {
        'long_name': 'Thermistor Temperature',
        'comment': ('Thermistor temperature refers to the internal instrument temperature of the pH sensor, as ' +
                    'measured by the thermistor. It is used to determine salinity and temperature dependent molar ' +
                    'absorptivities in the seawater sample in order to make an accurate pH estimation. This ' +
                    'variable represents the thermistor temperature measured at the end of the measurement cycle'),
        'units': 'degrees_Celsius',
        'ancillary_variables': 'raw_thermistor_end'
    },
    'seawater_ph': {
        'long_name': 'Seawater pH',
        'comment': ('pH is a measurement of the concentration of hydrogen ions in a solution. pH ranges from acidic ' +
                    'to basic on a scale from 0 to 14 with 7 being neutral.'),
        'standard_name': 'sea_water_ph_reported_on_total_scale',
        'data_product_identifier': 'PHWATER_L2',
        # 'units': ''    # deliberately left blank, no units for this value
        'ancillary_variables': ('blank_refrnc_434 blank_signal_434 blank_refrnc_578 blank_signal_578 ' +
                                'reference_434 signal_434 reference_578 signal_578 thermistor_temperature ' +
                                'practical_salinity')
    }
}


def phsen_datalogger(ds):
    """
    Takes PHSEN data recorded by the data loggers used in the CGSN/EA moorings and cleans up the data set to make
    it more user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that
    are of no use/value. Additionally, re-organize some of the variables to permit better assessments of the data.

    :param ds: initial PHSEN data set recorded by the data logger system and downloaded from OOI via the M2M system
    :return: cleaned up and reorganized data set
    """
    # drop some of the variables:
    #   passed_checksum == if it didn't, we wouldn't have a record, deleting
    #   record_type == there is only one, don't need this
    #   record_time == internal_timestamp, redundant so can remove
    #   dcl_controller_timestamp == time, redundant so can remove
    #   phsen_abcdef_signal_intensity_434, part of the light measurements array, redundant so can remove
    #   phsen_abcdef_signal_intensity_578, part of the light measurements array, redundant so can remove
    #   provenance == better to access with direct call to OOI M2M api, it doesn't work well in this format
    ds = ds.reset_coords()
    ds = ds.drop(['passed_checksum', 'record_type', 'record_time', 'phsen_abcdef_signal_intensity_434',
                  'phsen_abcdef_signal_intensity_578', 'dcl_controller_timestamp', 'provenance'])

    # convert the time values from a datetime64[ns] object to a floating point number with the time in seconds
    ds['internal_timestamp'] = ('time', dt64_epoch(ds.internal_timestamp))
    ds['internal_timestamp'].attrs = dict({
        'long_name': 'Internal SAMI-pH Clock Time',
        'standard_name': 'time',
        'units': 'seconds since 1970-01-01 00:00:00 0:00',
        'calendar': 'gregorian',
        'comment': ('Comparing the instrument internal clock versus the GPS referenced sampling time will allow for ' +
                    'calculations of the instrument clock offset and drift. Useful when working with the ' +
                    'recovered instrument data where no external GPS referenced clock is available.')
    })

    # rename some of the variables for better clarity
    rename = {
        'voltage_battery': 'raw_battery_voltage',
        'thermistor_start': 'raw_thermistor_start',
        'thermistor_end': 'raw_thermistor_end',
        'phsen_thermistor_temperature': 'thermistor_temperature',
        'phsen_abcdef_ph_seawater': 'seawater_ph',
        'phsen_abcdef_ph_seawater_qc_executed': 'seawater_ph_qc_executed',
        'phsen_abcdef_ph_seawater_qc_results': 'seawater_ph_qc_results'
    }
    ds = ds.rename(rename)

    # now we need to reset the light and reference arrays to named variables that will be more meaningful and useful in
    # the final data files
    nrec = len(ds['time'].values)
    light = np.array(np.vstack(ds['light_measurements'].values), dtype='int32')
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


def phsen_instrument(ds):
    """
    Takes PHSEN data recorded by internally by the instrument, and cleans up the data set to make it more
    user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that are
    of no use/value. Additionally, re-organize some of the variables to permit better assessments of the data.

    :param ds: initial PHSEN data set recorded internally by the instrument and downloaded from OOI via the M2M system
    :return: cleaned up data set
    """
    # drop some of the variables:
    #   record_type == there is only one, don't need this
    #   record_time == internal_timestamp == time, redundant so can remove both
    #   provenance == better to access with direct call to OOI M2M api, it doesn't work well in this format
    ds = ds.reset_coords()
    ds = ds.drop(['record_type', 'record_time', 'internal_timestamp', 'provenance'])

    # rename some of the variables for better clarity
    rename = {
        'voltage_battery': 'raw_battery_voltage',
        'thermistor_start': 'raw_thermistor_start',
        'thermistor_end': 'raw_thermistor_end',
        'phsen_thermistor_temperature': 'thermistor_temperature',
        'phsen_abcdef_ph_seawater': 'seawater_ph',
        'phsen_abcdef_ph_seawater_qc_executed': 'seawater_ph_qc_executed',
        'phsen_abcdef_ph_seawater_qc_results': 'seawater_ph_qc_results'
    }
    ds = ds.rename(rename)

    # now we need to reset the light and reference arrays to named variables that will be more meaningful and useful in
    # the final data files
    nrec = len(ds['time'].values)
    light = np.array(np.vstack(ds['light_measurements'].values), dtype='int32')
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
    data_types = ['deployment', 'raw_thermistor_end', 'raw_thermistor_start', 'raw_battery_voltage']
    for v in data_types:
        ds[v] = ds[v].astype('int32')

    return ds


def phsen_imodem(ds):
    """
    Takes PHSEN data recorded by over the inductive modem line, and cleans up the data set to make it more
    user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that are
    of no use/value. Additionally, re-organize some of the variables to permit better assessments of the data.

    :param ds: initial PHSEN data set recorded internally by the instrument and downloaded from OOI via the M2M system
    :return: cleaned up data set
    """
    # drop some of the variables:
    #   passed_checksum == useless, if it didn't we wouldn't have any data
    #   record_type == there is only one, don't need this
    #   record_time == internal_timestamp == time, redundant so can remove both
    #   phsen_abcdef_signal_intensity_434, part of the light measurements array, redundant so can remove
    #   phsen_abcdef_signal_intensity_578, part of the light measurements array, redundant so can remove
    #   provenance == better to access with direct call to OOI M2M api, it doesn't work well in this format
    ds = ds.reset_coords()
    ds = ds.drop(['passed_checksum', 'record_type', 'record_time', 'internal_timestamp', 'provenance',
                  'phsen_abcdef_signal_intensity_434', 'phsen_abcdef_signal_intensity_578'])

    # rename some of the variables for better clarity
    rename = {
        'voltage_battery': 'raw_battery_voltage',
        'thermistor_start': 'raw_thermistor_start',
        'thermistor_end': 'raw_thermistor_end',
        'phsen_thermistor_temperature': 'thermistor_temperature',
        'phsen_abcdef_ph_seawater': 'seawater_ph',
        'phsen_abcdef_ph_seawater_qc_executed': 'seawater_ph_qc_executed',
        'phsen_abcdef_ph_seawater_qc_results': 'seawater_ph_qc_results'
    }
    ds = ds.rename(rename)

    # now we need to reset the light and reference arrays to named variables that will be more meaningful and useful in
    # the final data files
    nrec = len(ds['time'].values)
    light = np.array(np.vstack(ds['light_measurements'].values), dtype='int32')
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
    data_types = ['deployment', 'raw_thermistor_end', 'raw_thermistor_start', 'raw_battery_voltage']
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

    # Request the data for download
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    if not r:
        exit_text = ('Request failed for %s-%s-%s. Check request.' % (site, node, sensor))
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
    if method in ['telemetered', 'recovered_host']:
        if re.match('.*imodem.*', stream):
            phsen = phsen_imodem(phsen)
        else:
            phsen = phsen_datalogger(phsen)
    else:
        phsen = phsen_instrument(phsen)

    vocab = get_vocabulary(site, node, sensor)[0]
    phsen = update_dataset(phsen, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(os.path.join(CONFIG['base_dir']['m2m_base'], args.outfile))
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    phsen.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
