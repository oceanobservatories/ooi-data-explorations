#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import numpy as np
import os

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, load_gc_thredds, \
    get_vocabulary, update_dataset, dict_update, ENCODINGS
from ooi_data_explorations.qartod.qc_processing import parse_qc

# load configuration settings
ATTRS = {
    'spectrum_average': {
        'long_name': 'Spectrum Average',
        'comment': 'Average value from the raw spectral measurements.',
        'units': 'counts'
    },
    'nitrate_concentration': {
        'long_name': 'Dissolved Nitrate Concentration',
        'comment': 'Dissolved Nitrate Concentration, uncorrected for temperature and salinity effects.',
        'units': 'umol L-1',
        'data_product_identifier': 'NITRDIS_L1',
    },
    'corrected_nitrate_concentration': {
        'long_name': 'Corrected Dissolved Nitrate Concentration',
        'standard_name': 'mole_concentration_of_nitrate_in_sea_water',
        'comment': 'Temperature and salinity corrected dissolved nitrate concentration.',
        'units': 'umol L-1',
        'data_product_identifier': 'NITRTSC_L2',
        'ancillary_variables': ('sea_water_temperature sea_water_practical_salinity raw_spectral_measurements '
                                'dark_value_used_for_fit')
    },
    'raw_spectral_measurements': {
        'long_name': 'Raw Spectral Measurements',
        'comment': ('The raw absorption spectra is an array of values, measured in counts and produced by the SUNA '
                    'based on the ultraviolet (UV) absorption characteristics of the components of seawater (including '
                    'dissolved nitrate). The raw spectral measurements, between 217 and 240 nm, are used to calculate '
                    'the nitrate concentrations.'),
        'data_product_identifier': 'NITROPT_L0',
        'units': 'counts'
    },
    'wavelength_index': {
        'long_name': 'Wavelength Index',
        'comment': ('Indexing array, with values between 0 and 255, for the raw spectral measurements. Each SUNA '
                    'measures between 190 and 370 nm, approximately every 0.7 nm. The exact wavelengths used are '
                    'specific to each instrument.'),
        # 'units': ''  # deliberately left blank, no units for this variable
    }
}


def suna_datalogger(ds, burst=True):
    """
    Takes SUNA data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some variables to permit better assessments of the data.

    :param ds: initial nutnr data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return ds: cleaned up data set
    """
    # drop some variables:
    #   checksum = used to parse data, is not parsed if the checksum fails so no longer needed
    #   frame_type = remove the dark frames if recorded, then remove
    #   humidity = not measured, no need to include
    ds = ds.reset_coords()
    # CGSN update: changed function order to handle arrays kept as bytes 
    ds['frame_type'] = ds['frame_type'].astype(str)
    ds = ds.where(ds.frame_type == 'SLF', drop=True)  # remove the dark frames
    ds = ds.drop(['checksum', 'frame_type', 'humidity'])

    # convert the instrument date and time values to a floating point number with the time in seconds, and then
    # update the internal_timestamp with the values as it was not set correctly during parsing and defaults to 1900.
    date_of_sample = ds.date_of_sample.values.astype(str)
    ddate = [datetime.datetime.strptime(dd[:7], '%Y%j') for dd in date_of_sample]
    dhour = [datetime.timedelta(hours=dh) for dh in ds.time_of_sample.values]
    internal_timestamp = np.array([np.datetime64(ddate[i] + dhour[i]) for i in range(len(ddate))])
    ds['internal_timestamp'] = ('time', internal_timestamp.astype(float) / 10.0 ** 6)
    ds['internal_timestamp'].attrs = dict({
        'long_name': 'Internal NUTNR Clock Time',
        'standard_name': 'time',
        'units': 'seconds since 1900-01-01T00:00:00.000Z',
        'calendar': 'gregorian',
        'comment': ('Comparing the instrument internal clock versus the GPS referenced sampling time will allow for '
                    'calculations of the instrument clock offset and drift. Useful when working with the '
                    'recovered instrument data where no external GPS referenced clock is available.')
    })
    ds = ds.drop(['date_of_sample', 'time_of_sample'])

    # check for data from a co-located CTD, if not present add with appropriate attributes
    if 'sea_water_temperature' not in ds.variables:
        ds['sea_water_temperature'] = ('time', ds['deployment'].data * np.nan)
        ds['sea_water_temperature'].attrs = {
            'comment': ('Normally this would be sea water temperature data from a co-located CTD. However, data from ' +
                        'that sensor is unavailable. This value has been filled with NaNs to preserve the structure ' +
                        'of the data set.'),
            'data_product_identifier': 'TEMPWAT_L1',
            'long_name': 'Sea Water Temperature',
            'standard_name': 'sea_water_temperature',
            'units': 'degree_Celsius'
        }

        ds['sea_water_practical_salinity'] = ('time', ds['deployment'].data * np.nan)
        ds['sea_water_practical_salinity'].attrs = {
            'long_name': 'Sea Water Practical Salinity',
            'standard_name': 'sea_water_practical_salinity',
            'units': '1',
            'comment': ('Normally this would be seawater salinity data from a co-located CTD. However, data from ' +
                        'that sensor is unavailable. This value has been filled with NaNs to preserve the structure ' +
                        'of the data set.'),
            'data_product_identifier': 'PRACSAL_L2'
        }

    # rename some variables for better clarity
    rename = {
        'nutnr_absorbance_at_254_nm': 'absorbance_at_254_nm',
        'nutnr_absorbance_at_350_nm': 'absorbance_at_350_nm',
        'nutnr_bromide_trace': 'bromide_trace',
        'nutnr_current_main': 'current_main',
        'nutnr_dark_value_used_for_fit': 'dark_value_used_for_fit',
        'nutnr_fit_base_1': 'fit_base_1',
        'nutnr_fit_base_2': 'fit_base_2',
        'nutnr_fit_rmse': 'fit_rmse',
        'nutnr_integration_time_factor': 'integration_time_factor',
        'nutnr_nitrogen_in_nitrate': 'nitrogen_in_nitrate',
        'nutnr_spectrum_average': 'spectrum_average',
        'nutnr_voltage_int': 'voltage_instrument',
        'temp_spectrometer': 'temperature_spectrometer',
        'temp_lamp': 'temperature_lamp',
        'temp_interior': 'temperature_interior',
        'spectral_channels': 'raw_spectral_measurements',
        'salinity_corrected_nitrate': 'corrected_nitrate_concentration',
        'salinity_corrected_nitrate_qc_results': 'corrected_nitrate_concentration_qc_results',
        'salinity_corrected_nitrate_qc_executed': 'corrected_nitrate_concentration_qc_executed',
        'wavelength': 'wavelength_index'
    }
    for key, value in rename.items():
        if key in ds.variables:
            ds = ds.rename({key: value})
            ds[value].attrs['ooinet_variable_name'] = key

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # address incorrectly set units and variable types
    ds['dark_value_used_for_fit'].attrs['units'] = 'counts'
    # CGSN update: switched to simple type conversions
    ds['serial_number'] = ds['serial_number'].astype(int)
    #ds['serial_number'].attrs = dict({
    #    'long_name': 'Serial Number',
        # 'units': '', deliberately left blank, unitless value
    #    'comment': ('Instrument serial number'),
    #})

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    if burst:   # re-sample the data to a defined time interval using a median average
        # create the burst averaging
        burst = ds
        burst.load()
        burst = burst.resample(time='900s', base=3150, loffset='450s', skipna=True).median(keep_attrs=True)
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # save the newly averaged data
        ds = burst

        # and reset some data types
        data_types = ['deployment', 'spectrum_average', 'serial_number', 'dark_value_used_for_fit',
                      'raw_spectral_measurements']
        for v in data_types:
            ds[v] = ds[v].astype('int32')

    return ds


def suna_instrument(ds, burst=True):
    """
    Takes SUNA data recorded internally by the Sea-Bird SUNA dissolved nitrate
    sensor and cleans up the data set to make it more user-friendly.  Primary
    task is renaming parameters and dropping some that are of limited use.
    Additionally, re-organize some variables to permit better assessments of
    the data.

    :param ds: initial nutnr data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return ds: cleaned up data set
    """
    # drop some variables:
    #   checksum = used to parse data, is not parsed if the checksum fails so no longer needed
    #   frame_type = remove the dark frames if recorded, then remove
    #   humidity = not measured, no need to include
    #   internal_timestamp == time, redundant so can remove
    #   date_of_sample = used to construct the internal_timestamp
    #   time_of_sample = used to construct the internal_timestamp
    ds = ds.reset_coords()
    # CGSN update: switched to simple type conversion - handles byte type arrays
    ds['frame_type'] = ('time', [''.join(x.astype(str)) for x in ds.frame_type.data])
    ds = ds.where(ds.frame_type == 'SLF', drop=True)  # remove the dark frames
    ds = ds.drop(['checksum', 'frame_type', 'humidity', 'date_of_sample', 'time_of_sample'])

    # check for data from a co-located CTD, if not present add with appropriate attributes
    if 'sea_water_temperature' not in ds.variables:
        ds['sea_water_temperature'] = ('time', ds['deployment'].data * np.nan)
        ds['sea_water_temperature'].attrs = {
            'comment': ('Normally this would be sea water temperature data from a co-located CTD. However, data from ' +
                        'that sensor is unavailable. This value has been filled with NaNs to preserve the structure ' +
                        'of the data set.'),
            'data_product_identifier': 'TEMPWAT_L1',
            'long_name': 'Sea Water Temperature',
            'standard_name': 'sea_water_temperature',
            'units': 'degree_Celsius'
        }

        ds['sea_water_practical_salinity'] = ('time', ds['deployment'].data * np.nan)
        ds['sea_water_practical_salinity'].attrs = {
            'long_name': 'Sea Water Practical Salinity',
            'standard_name': 'sea_water_practical_salinity',
            'units': '1',
            'comment': ('Normally this would be seawater salinity data from a co-located CTD. However, data from ' +
                        'that sensor is unavailable. This value has been filled with NaNs to preserve the structure ' +
                        'of the data set.'),
            'data_product_identifier': 'PRACSAL_L2'
        }

    # rename some variables for better clarity
    rename = {
        'nutnr_absorbance_at_254_nm': 'absorbance_at_254_nm',
        'nutnr_absorbance_at_350_nm': 'absorbance_at_350_nm',
        'nutnr_bromide_trace': 'bromide_trace',
        'nutnr_current_main': 'current_main',
        'nutnr_dark_value_used_for_fit': 'dark_value_used_for_fit',
        'nutnr_fit_base_1': 'fit_base_1',
        'nutnr_fit_base_2': 'fit_base_2',
        'nutnr_fit_rmse': 'fit_rmse',
        'nutnr_integration_time_factor': 'integration_time_factor',
        'nutnr_nitrogen_in_nitrate': 'nitrogen_in_nitrate',
        'nutnr_spectrum_average': 'spectrum_average',
        'nutnr_voltage_int': 'voltage_instrument',
        'temp_spectrometer': 'temperature_spectrometer',
        'temp_lamp': 'temperature_lamp',
        'temp_interior': 'temperature_interior',
        'spectral_channels': 'raw_spectral_measurements',
        'salinity_corrected_nitrate': 'corrected_nitrate_concentration',
        'salinity_corrected_nitrate_qc_results': 'corrected_nitrate_concentration_qc_results',
        'salinity_corrected_nitrate_qc_executed': 'corrected_nitrate_concentration_qc_executed',
        'wavelength': 'wavelength_index'
    }
    for key, value in rename.items():
        if key in ds.variables:
            ds = ds.rename({key: value})
            ds[value].attrs['ooinet_variable_name'] = key

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # address incorrectly set units and variable types
    ds['dark_value_used_for_fit'].attrs['units'] = 'counts'
    # CGSN update: switched to type conversions
    ds['serial_number'] = ds['serial_number'].astype(int)
    #ds['serial_number'].attrs = dict({
    #    'long_name': 'Serial Number',
        # 'units': '', deliberately left blank, unitless value
    #    'comment': ('Instrument serial number'),
    #})

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    if burst:   # re-sample the data to a defined time interval using a median average
        # create the burst averaging
        burst = ds
        burst.load()
        burst = burst.resample(time='900s', base=3150, loffset='450s', skipna=True).median(keep_attrs=True)
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # save the newly averaged data
        ds = burst

        # and reset some data types
        data_types = ['deployment', 'spectrum_average', 'serial_number', 'dark_value_used_for_fit',
                      'raw_spectral_measurements']
        for v in data_types:
            ds[v] = ds[v].astype('int32')

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

    if stream not in ['suna_dcl_recovered', 'suna_instrument_recovered']:
        raise SystemExit('Data from the Satlantic ISUS dissolved nitrate sensor will not be supported at this time.')

    # check if we are specifying a deployment or a specific date and time range
    if not deploy or (start and stop):
        return SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')

    # if we are specifying a deployment number, then get the data from the Gold Copy THREDDS server
    if deploy:
        # download the data for the deployment
        nutnr = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*NUTNR.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not nutnr:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, deployment %d.' % (site, node, sensor, method,
                                                                                    stream, deploy))
            raise SystemExit(exit_text)
    else:
        # otherwise, request the data for download from OOINet via the M2M API using the specified dates
        r = m2m_request(site, node, sensor, method, stream, start, stop)
        if not r:
            exit_text = ('Request failed for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                  stream, start, stop))
            raise SystemExit(exit_text)

        # Valid M2M request, start downloading the data
        nutnr = m2m_collect(r, '.*NUTNR.*\\.nc$')

        # check to see if we downloaded any data
        if not nutnr:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize
    if method in ['telemetered', 'recovered_host']:
        nutnr = suna_datalogger(nutnr, burst)
    else:
        nutnr = suna_instrument(nutnr, burst)

    vocab = get_vocabulary(site, node, sensor)[0]
    nutnr = update_dataset(nutnr, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    nutnr.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
