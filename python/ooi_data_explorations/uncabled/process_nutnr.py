#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import numpy as np
import os

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, get_deployment_dates, \
    get_vocabulary, update_dataset, dict_update, ENCODINGS

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
        'long_name': 'Corrected Nitrate Concentration',
        'standard_name': 'mole_concentration_of_nitrate_in_sea_water',
        'comment': 'Temperature and salinity corrected dissolved nitrate concentration.',
        'units': 'umol L-1',
        'data_product_identifier': 'NITRTSC_L2',
        'ancillary_variables': ('seawater_temperature practical_salinity raw_spectral_measurements '
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


def nutnr_datalogger(ds, burst=True):
    """
    Takes nutnr data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial nutnr data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   checksum = used to parse data, is not parsed if the checksum fails so no longer needed
    #   frame_type = only frame types are SLF per the parser, don't need to repeat
    #   humidity = not measured, no need to include
    ds = ds.reset_coords()
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
        'units': 'seconds since 1970-01-01 00:00:00 0:00',
        'calendar': 'gregorian',
        'comment': ('Comparing the instrument internal clock versus the GPS referenced sampling time will allow for '
                    'calculations of the instrument clock offset and drift. Useful when working with the '
                    'recovered instrument data where no external GPS referenced clock is available.')
    })
    ds = ds.drop(['date_of_sample', 'time_of_sample'])

    # check for data from a co-located CTD, if not present add with appropriate attributes
    if 'temp' not in ds.variables:
        ds['temp'] = ('time', ds['deployment'] * np.nan)
        ds['temp'].attrs = {
            'comment': ('Normally this would be seawater temperature data from a co-located CTD. However, data from ' +
                        'that sensor is unavailable. This value has been filled with NaNs to preserve the structure ' +
                        'of the data set.'),
            'data_product_identifier': 'TEMPWAT_L1',
            'long_name': 'Seawater Temperature',
            'standard_name': 'sea_water_temperature',
            'units': 'degree_Celsius'
        }

        ds['practical_salinity'] = ('time', ds['deployment'] * np.nan)
        ds['practical_salinity'].attrs = {
            'long_name': 'Practical Salinity',
            'standard_name': 'sea_water_practical_salinity',
            'units': '1',
            'comment': ('Normally this would be seawater salinity data from a co-located CTD. However, data from ' +
                        'that sensor is unavailable. This value has been filled with NaNs to preserve the structure ' +
                        'of the data set.'),
            'data_product_identifier': 'PRACSAL_L2'
        }

    # rename some of the variables for better clarity
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
        'temp': 'seawater_temperature',
        'spectral_channels': 'raw_spectral_measurements',
        'salinity_corrected_nitrate': 'corrected_nitrate_concentration',
        'salinity_corrected_nitrate_qc_results': 'corrected_nitrate_concentration_qc_results',
        'salinity_corrected_nitrate_qc_executed': 'corrected_nitrate_concentration_qc_executed',
        'wavelength': 'wavelength_index'
    }
    ds = ds.rename(rename)

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # add the original variable name as an attribute, if renamed
    for key, value in rename.items():
        ds[value].attrs['ooinet_variable_name'] = key

    # address incorrectly set units
    ds['seawater_temperature'].attrs['units'] = 'degree_Celsius'
    ds['temperature_lamp'].attrs['units'] = 'degree_Celsius'
    ds['temperature_interior'].attrs['units'] = 'degree_Celsius'
    ds['temperature_spectrometer'].attrs['units'] = 'degree_Celsius'
    ds['nitrate_concentration'].attrs['units'] = 'umol L-1'
    ds['corrected_nitrate_concentration'].attrs['units'] = 'umol L-1'
    ds['dark_value_used_for_fit'].attrs['units'] = 'counts'
    ds['serial_number'] = ds['serial_number'].astype('int32')

    if burst:   # re-sample the data to a defined time interval using a median average
        # create the burst averaging
        burst = ds
        burst = burst.resample(time='900s', base=3150, loffset='450s', keep_attrs=True, skipna=True).median()
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # reset the attributes...which keep_attrs should do...
        burst.attrs = ds.attrs
        for v in burst.variables:
            burst[v].attrs = ds[v].attrs

        # save the newly average data
        ds = burst

    # and reset some of the data types
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

    if stream not in ['suna_dcl_recovered']:
        exit_text = ('Currently the only stream supported is suna_dcl_recovered, you requested %s.' % stream)
        raise SystemExit(exit_text)

    # Request the data for download
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    if not r:
        exit_text = ('Request failed for %s-%s-%s. Check request.' % (site, node, sensor))
        raise SystemExit(exit_text)

    # Valid request, start downloading the data
    if deploy:
        nutnr = m2m_collect(r, ('.*deployment%04d.*NUTNR.*\\.nc$' % deploy))
    else:
        nutnr = m2m_collect(r, '.*NUTNR.*\\.nc$')

    if not nutnr:
        exit_text = ('Data unavailable for %s-%s-%s. Check request.' % (site, node, sensor))
        raise SystemExit(exit_text)

    # clean-up and reorganize
    nutnr = nutnr_datalogger(nutnr, burst)
    vocab = get_vocabulary(site, node, sensor)[0]
    nutnr = update_dataset(nutnr, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    nutnr.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
