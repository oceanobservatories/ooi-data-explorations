#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
import os
import pandas as pd
import xarray as xr

from scipy.interpolate import griddata

from ooi_data_explorations.common import inputs, load_gc_thredds, m2m_collect, m2m_request, get_vocabulary, \
    update_dataset, ENCODINGS
from ooi_data_explorations.qartod.qc_processing import parse_qc

# load configuration settings
FILL_INT = -9999999
ATTRS = dict({
    # dataset attributes --> raw values
    'raw_irradiance_412': {
        'long_name': 'Raw Downwelling Irradiance 412 nm',
        'comment': ('Raw downwelling spectral irradiance at 412 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-412',
    },
    'raw_irradiance_444': {
        'long_name': 'Raw Downwelling Irradiance 444 nm',
        'comment': ('Raw downwelling spectral irradiance at 444 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-444',
    },
    'raw_irradiance_490': {
        'long_name': 'Raw Downwelling Irradiance 490 nm',
        'comment': ('Raw downwelling spectral irradiance at 490 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-490',
    },
    'raw_irradiance_510': {
        'long_name': 'Raw Downwelling Irradiance 510 nm',
        'comment': ('Raw downwelling spectral irradiance at 510 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-510',
    },
    'raw_irradiance_555': {
        'long_name': 'Raw Downwelling Irradiance 555 nm',
        'comment': ('Raw downwelling spectral irradiance at 555 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-555',
    },
    'raw_irradiance_620': {
        'long_name': 'Raw Downwelling Irradiance 620 nm',
        'comment': ('Raw downwelling spectral irradiance at 620 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-620',
    },
    'raw_irradiance_683': {
        'long_name': 'Raw Downwelling Irradiance 683 nm',
        'comment': ('Raw downwelling spectral irradiance at 683 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-683',
    },
    # dataset attributes --> derived values
    'downwelling_irradiance_412': {
        'long_name': 'Downwelling Spectral Irradiance at 412 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 412,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 412 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-412',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_412'
    },
    'downwelling_irradiance_444': {
        'long_name': 'Downwelling Spectral Irradiance at 444 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 444,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 444 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-444',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_444'
    },
    'downwelling_irradiance_490': {
        'long_name': 'Downwelling Spectral Irradiance at 490 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 490,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 490 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-490',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_490'
    },
    'downwelling_irradiance_510': {
        'long_name': 'Downwelling Spectral Irradiance at 510 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 510,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 510 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-510',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_510'
    },
    'downwelling_irradiance_555': {
        'long_name': 'Downwelling Spectral Irradiance at 555 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 555,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 555 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-555',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_555'
    },
    'downwelling_irradiance_620': {
        'long_name': 'Downwelling Spectral Irradiance at 620 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 620,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 620 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-620',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_620'
    },
    'downwelling_irradiance_683': {
        'long_name': 'Downwelling Spectral Irradiance at 683 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 683,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 683 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-683',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_683'
    }
})


def quality_checks(ds):
    """
    Assessment of the raw data and the calculated parameters for quality
    using a susbset of the QARTOD flags to indicate the quality. QARTOD
    flags used are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail

    The final flag value represents the worst case assessment of the data
    quality.

    :param ds: xarray dataset with the raw signal data and the calculated
               bio-optical parameters
    :return beta_flag: QARTOD quality flags for the backscatter measurements
    :return cdom_flag: QARTOD quality flags for the CDOM measurements
    :return chl_flag: QARTOD quality flags for the chlorophyll measurements
    """
    return None


def spkir_datalogger(ds, burst=False):
    """
    Takes spkir data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly. Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial spkir data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   instrument_id == captured in global metadata (always the OCR-507)
    #   dcl_controller_timestamp == empty string
    #   passed_checksum == not ingested if it doesn't pass, don't need to restate that here
    #   internal_timestamp == superseded by time, redundant so can remove
    #   sample_delay == superseded by time, redundant so can remove
    #   timer == superseded by time, redundant so can remove
    drop_vars = ['instrument_id', 'dcl_controller_timestamp', 'passed_checksum',
                 'internal_timestamp', 'sample_delay', 'timer']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # rename some parameters here to get a better defined data set with cleaner attributes
    rename = {
        'vin_sense': 'input_voltage',
        'va_sense': 'analog_rail_voltage'
    }
    ds = ds.rename(rename)

    # convert raw voltages and the instrument temperature to engineering units
    ds['input_voltage'] = ds['input_voltage'] * 0.03
    ds['analog_rail_voltage'] = ds['analog_rail_voltage'] * 0.03
    ds['internal_temperature'] = -50 + ds['internal_temperature'] * 0.5

    # add the raw downwelling spectral irradiance measurements to the data set split out by wavelength
    raw = ds['channel_array']
    ds['raw_irradiance_412'] = raw[:, 0]
    ds['raw_irradiance_444'] = raw[:, 1]
    ds['raw_irradiance_490'] = raw[:, 2]
    ds['raw_irradiance_510'] = raw[:, 3]
    ds['raw_irradiance_555'] = raw[:, 4]
    ds['raw_irradiance_620'] = raw[:, 5]
    ds['raw_irradiance_683'] = raw[:, 6]

    # add the converted downwelling spectral irradiance measurements to the data set split out by wavelength
    ed = ds['spkir_abj_cspp_downwelling_vector']
    ds['downwelling_irradiance_412'] = ed[:, 0]
    ds['downwelling_irradiance_444'] = ed[:, 1]
    ds['downwelling_irradiance_490'] = ed[:, 2]
    ds['downwelling_irradiance_510'] = ed[:, 3]
    ds['downwelling_irradiance_555'] = ed[:, 4]
    ds['downwelling_irradiance_620'] = ed[:, 5]
    ds['downwelling_irradiance_683'] = ed[:, 6]

    # remove the original 2D variables
    ds = ds.drop_vars(['spectra', 'channel_array', 'spkir_abj_cspp_downwelling_vector'])

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    if burst:
        # re-sample the data collected in burst mode using a 15-minute median average
        burst = ds.resample(time='900s', skipna=True).median(dim='time', keep_attrs=True)

        # remove the missing rows
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # save the newly average data
        ds = burst

    return ds


def spkir_cspp(ds):
    """
    Takes SPKIR data recorded by the CSPP loggers used by the Endurance Array
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial SPKIR data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   seawater_scattering_coefficient == not used
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl'])

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'pressure': 'seawater_pressure',
        'pressure_qc_executed': 'seawater_pressure_qc_executed',
        'pressure_qc_results': 'seawater_pressure_qc_results',
        'temperature': 'seawater_temperature',
        'salinity': 'practical_salinity',
        'raw_signal_chl': 'raw_chlorophyll',
        'fluorometric_chlorophyll_a': 'estimated_chlorophyll',
        'fluorometric_chlorophyll_a_qc_executed': 'estimated_chlorophyll_qc_executed',
        'fluorometric_chlorophyll_a_qc_results': 'estimated_chlorophyll_qc_results',
        'raw_signal_cdom': 'raw_cdom',
        'raw_signal_beta': 'raw_backscatter',
        'total_volume_scattering_coefficient': 'beta_700',
        'total_volume_scattering_coefficient_qc_executed': 'beta_700_qc_executed',
        'total_volume_scattering_coefficient_qc_results': 'beta_700_qc_results',
        'optical_backscatter': 'bback',
        'optical_backscatter_qc_executed': 'bback_qc_executed',
        'optical_backscatter_qc_results': 'bback_qc_results',
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

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    # create qc flags for the data and add them to the OOI qc flags
    beta_flag, cdom_flag, chl_flag = quality_checks(ds)
    ds['beta_700_qc_summary_flag'] = ('time', (np.array([ds.beta_700_qc_summary_flag,
                                                         beta_flag])).max(axis=0, initial=1))
    ds['fluorometric_cdom_qc_summary_flag'] = ('time', (np.array([ds.fluorometric_cdom_qc_summary_flag,
                                                                 cdom_flag])).max(axis=0, initial=1))
    ds['estimated_chlorophyll_qc_summary_flag'] = ('time', (np.array([ds.estimated_chlorophyll_qc_summary_flag,
                                                                      chl_flag])).max(axis=0, initial=1))

    return ds


def spkir_wfp(ds, grid=False):
    """
    Takes SPKIR data recorded by the Wire-Following Profilers (used by CGSN/EA
    as part of the coastal and global arrays) and cleans up the data set to
    make it more user-friendly.  Primary task is renaming parameters and
    dropping some that are of limited use. Additionally, re-organize some of
    the variables to permit better assessments of the data.

    :param ds: initial SPKIR data set downloaded from OOI via the M2M system
    :param grid: boolean flag for whether the data should be gridded
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   seawater_scattering_coefficient == not used
    #   raw_internal_temp == not available, NaN filled
    ds = ds.reset_coords()
    drop_vars = ['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                 'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'raw_internal_temp']
    for v in drop_vars:
        if v in ds.variables:
            ds = ds.drop_vars(v)

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'int_ctd_pressure': 'seawater_pressure',
        'ctdpf_ckl_seawater_temperature': 'seawater_temperature',
        'raw_signal_chl': 'raw_chlorophyll',
        'fluorometric_chlorophyll_a': 'estimated_chlorophyll',
        'fluorometric_chlorophyll_a_qc_executed': 'estimated_chlorophyll_qc_executed',
        'fluorometric_chlorophyll_a_qc_results': 'estimated_chlorophyll_qc_results',
        'raw_signal_cdom': 'raw_cdom',
        'raw_signal_beta': 'raw_backscatter',
        'total_volume_scattering_coefficient': 'beta_700',
        'total_volume_scattering_coefficient_qc_executed': 'beta_700_qc_executed',
        'total_volume_scattering_coefficient_qc_results': 'beta_700_qc_results',
        'optical_backscatter': 'bback',
        'optical_backscatter_qc_executed': 'bback_qc_executed',
        'optical_backscatter_qc_results': 'bback_qc_results',
    }
    for key in rename.keys():
        if key in ds.variables:
            ds = ds.rename({key: rename.get(key)})

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # add the original variable name as an attribute, if renamed
    for key, value in rename.items():
        ds[value].attrs['ooinet_variable_name'] = key

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    # create qc flags for the data and add them to the OOI qc flags
    beta_flag, cdom_flag, chl_flag = quality_checks(ds)
    ds['beta_700_qc_summary_flag'] = ('time', (np.array([ds.beta_700_qc_summary_flag,
                                                         beta_flag])).max(axis=0, initial=1))
    ds['fluorometric_cdom_qc_summary_flag'] = ('time', (np.array([ds.fluorometric_cdom_qc_summary_flag,
                                                                 cdom_flag])).max(axis=0, initial=1))
    ds['estimated_chlorophyll_qc_summary_flag'] = ('time', (np.array([ds.estimated_chlorophyll_qc_summary_flag,
                                                                      chl_flag])).max(axis=0, initial=1))

    if grid:
        # clear out any duplicate time stamps
        _, index = np.unique(ds['time'], return_index=True)
        ds = ds.isel(time=index)

        # since the scipy griddata function cannot use the time values as is (get converted to nanoseconds, which
        # is too large of a value), we need to temporarily convert them to a floating point number in days since
        # the start of the data record; we can then use that temporary date/time array for the gridding.
        base_time = ds['time'].min().values
        dt = (ds['time'] - base_time).astype(float) / 1e9 / 60 / 60 / 24

        # construct the new grid, using 1 m depth bins from 30 to 510 m, and daily intervals from the start of
        # the record to the end (centered on noon UTC).
        depth_range = np.arange(30, 511, 1)
        time_range = np.arange(0.5, np.ceil(dt.max()) + 0.5, 1)
        gridded_time = base_time.astype('M8[D]') + pd.to_timedelta(time_range, unit='D')

        # grid the data, adding the results to a list of data arrays
        gridded = []
        for v in ds.variables:
            if v not in ['time', 'depth']:
                # grid the data for each variable
                gdata = griddata((dt.values, ds['depth'].values), ds[v].values,
                                 (time_range[None, :], depth_range[:, None]),
                                 method='linear')

                # add the data to a data array
                da = xr.DataArray(name=v, data=gdata, coords=[("depth", depth_range), ("time", gridded_time)])
                da.attrs = ds[v].attrs

                # reset the data types and fill values for floats and ints
                if ds[v].dtype == np.dtype(int):
                    da = da.where(np.isnan is True, FILL_INT)
                    da.attrs['_FillValue'] = FILL_INT
                    da = da.astype(int)
                else:
                    da.attrs['_FillValue'] = np.nan
                    da = da.astype(float)

                # add to the list
                gridded.append(da)

        # recombine the gridded data arrays into a single dataset
        gridded = xr.merge(gridded)
        gridded.attrs = ds.attrs
        ds = gridded

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

    # check if we are specifying a deployment or a specific date and time range
    if not deploy or (start and stop):
        return SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')

    # if we are specifying a deployment number, then get the data from the Gold Copy THREDDS server
    if deploy:
        # download the data for the deployment
        spkir = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*SPKIR.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not spkir:
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
        spkir = m2m_collect(r, '.*SPKIR.*\\.nc$')

        # check to see if we downloaded any data
        if not spkir:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize the data
    if node == 'SP001':
        # this SPKIR is part of a CSPP
        spkir = spkir_cspp(spkir)
    elif node == 'WFP01':
        # this SPKIR is part of a Wire-Following Profiler
        spkir = spkir_wfp(spkir)
    elif node == 'SBD17':
        # this SPKIR is connected to the CTDBP on an EA Inshore Surface Mooring
        spkir = spkir_instrument(spkir)
        if not spkir:
            # there was no data after removing all the 0's
            sys.exit()
    else:
        # this SPKIR is stand-alone on one of the moorings
        spkir = spkir_datalogger(spkir, burst)

    vocab = get_vocabulary(site, node, sensor)[0]
    spkir = update_dataset(spkir, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    spkir.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
