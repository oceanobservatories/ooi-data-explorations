#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
import os
import pandas as pd
import warnings
import xarray as xr

from scipy.interpolate import griddata

from ooi_data_explorations.common import inputs, load_gc_thredds, m2m_collect, m2m_request, get_vocabulary, \
    update_dataset, FILL_INT, FILL_FLOAT, ENCODINGS
from ooi_data_explorations.profilers import create_profile_id
from ooi_data_explorations.qartod.qc_processing import parse_qc

# load configuration settings
ATTRS = dict({
    'raw_backscatter': {
        'long_name': 'Raw Optical Backscatter at 700 nm',
        'units': 'counts',
        'comment': 'Raw optical backscatter measurements at 700 nm.',
        'data_product_identifier': 'FLUBSCT_L0'
    },
    'raw_chlorophyll': {
        'long_name': 'Raw Chlorophyll Fluorescence',
        'units': 'counts',
        'comment': 'Raw chlorophyll fluorescence (470 nm excitation/695 nm emission) measurements.',
        'data_product_identifier': 'CHLAFLO_L0'
    },
    'raw_fluorometric_cdom': {
        'long_name': 'Raw CDOM Fluorescence',
        'units': 'counts',
        'comment': 'Raw CDOM fluorescence (370 nm excitation/460 nm emission) measurements.',
        'data_product_identifier': 'CDOMFLO_L0'
    },
    'estimated_chlorophyll': {
        'long_name': 'Estimated Chlorophyll Concentration',
        'standard_name': 'mass_concentration_of_chlorophyll_in_sea_water',
        'units': 'ug L-1',
        'comment': ('Estimated chlorophyll concentration based upon a calibration curve derived from a fluorescent '
                    'proxy approximately equal to 25 ug/l of a Thalassiosira weissflogii phytoplankton culture. This '
                    'measurement is considered to be an estimate only of the true chlorophyll concentration.'),
        'data_product_identifier': 'CHLAFLO_L1',
        'ancillary_variables': 'raw_chlorophyll estimated_chlorophyll_qc_executed estimated_chlorophyll_qc_results'
    },
    'fluorometric_cdom': {
        'long_name': 'Fluorometric CDOM Concentration',
        'standard_name': ('concentration_of_colored_dissolved_organic_matter_in_sea_water_expressed_as_equivalent'
                          '_mass_fraction_of_quinine_sulfate_dihydrate'),
        'units': 'ppb',
        'comment': ('More commonly referred to as Chromophoric Dissolved Organic Matter (CDOM). CDOM plays an '
                    'important role in the carbon cycling and biogeochemistry of coastal waters. It occurs '
                    'naturally in aquatic environments primarily as a result of tannins released from decaying '
                    'plant and animal matter, and can enter coastal areas in river run-off containing organic '
                    'materials leached from soils.'),
        'data_product_identifier': 'CDOMFLO_L1',
        'ancillary_variables': 'raw_fluorometric_cdom fluorometric_cdom_qc_executed fluorometric_cdom_qc_results'
    },
    'beta_700': {
        'long_name': 'Volume Scattering Function at 700 nm',
        'standard_name': 'volume_scattering_function_of_radiative_flux_in_sea_water',
        'units': 'm-1 sr-1',
        'radiation_wavelength': 700.0,
        'radiation_wavelength_unit': 'nm',
        'comment': ('Radiative flux is the sum of shortwave and longwave radiative fluxes. Scattering of '
                    'radiation is its deflection from its incident path without loss of energy. The volume '
                    'scattering function is the intensity (flux per unit solid angle) of scattered radiation per '
                    'unit length of scattering medium, normalised by the incident radiation flux.'),
        'data_product_identifier': 'FLUBSCT_L1',
        'ancillary_variables': 'raw_backscatter beta_700_qc_executed beta_700_qc_results'
    },
    'bback': {
        'long_name': 'Total Optical Backscatter at 700 nm',
        'standard_name': 'volume_backwards_scattering_coefficient_of_radiative_flux_in_sea_water',
        'units': 'm-1',
        'comment': ('Total (particulate + water) optical backscatter at 700 nm, derived from the Volume '
                    'Scattering Function and corrected for effects of temperature and salinity.'),
        'data_product_identifier': 'FLUBSCT_L2',
        'ancillary_variables': ('beta_700 sea_water_temperature sea_water_practical_salinity '
                                'sea_water_scattering_coefficient bback_qc_executed bback_qc_results'),
        '_FillValue': FILL_FLOAT
    },
    'sea_water_scattering_coefficient': {
        'long_name': 'Sea Water Optical Backscatter at 700 nm',
        'units': 'm-1',
        'comment': ('Theoretical estimation of the optical backscatter for pure seawater at 700 nm adjusted for the'
                    'effects of temperature and salinity. This value is added to the particulate optical backscatter '
                    'measurement to create the total optical backscatter measurement contained in this data set.'),
        'ancillary_variables': 'sea_water_temperature sea_water_practical_salinity',
        '_FillValue': FILL_FLOAT
    },
    'sea_water_practical_salinity': {
        'long_name': 'Practical Salinity',
        'standard_name': 'sea_water_practical_salinity',
        'units': '1',
        'comment': ('Salinity is generally defined as the concentration of dissolved salt in a parcel of sea water. ' 
                    'Practical Salinity is a more specific unitless quantity calculated from the conductivity of ' 
                    'sea water and adjusted for temperature and pressure. It is approximately equivalent to Absolute ' 
                    'Salinity (the mass fraction of dissolved salt in sea water), but they are not interchangeable. '
                    'Measurements are from a co-located CTD.'),
        'data_product_identifier': 'PRACSAL_L2',
        '_FillValue': FILL_FLOAT
    },
    'sea_water_temperature': {
        'long_name': 'Sea Water Temperature',
        'standard_name': 'sea_water_temperature',
        'units': 'degree_Celsius',
        'comment': ('Sea water temperature is the in situ temperature of the sea water. Measurements are from a '
                    'co-located CTD'),
        'data_product_identifier': 'TEMPWAT_L1',
        '_FillValue': FILL_FLOAT
    },
})


def quality_checks(ds):
    """
    Assessment of the raw data and the calculated parameters for quality
    using a subset of the QARTOD flags to indicate the quality. QARTOD
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
    max_counts = 4115   # counts should be greater than 0 and less than 4120 +/- 5
    beta_flag = ds['time'].astype('int32') * 0 + 1   # default flag values, no errors
    cdom_flag = ds['time'].astype('int32') * 0 + 1   # default flag values, no errors
    chl_flag = ds['time'].astype('int32') * 0 + 1   # default flag values, no errors

    # test the min/max values of the raw measurements
    m = (ds.raw_backscatter <= 0) | (ds.raw_backscatter > max_counts)
    beta_flag[m] = 4    # raw volume scattering coefficient values off scale
    m = (ds.raw_fluorometric_cdom <= 0) | (ds.raw_fluorometric_cdom > max_counts)
    cdom_flag[m] = 4    # raw CDOM values off scale
    m = (ds.raw_chlorophyll <= 0) | (ds.raw_chlorophyll > max_counts)
    chl_flag[m] = 4     # raw chlorophyll values off scale

    # test the min/max values of the derived measurements (values from the vendor documentation)
    m = (ds.bback <= 0) | (ds.bback > 3)  # scattering measurement range
    beta_flag[m] = 4
    m = (ds.fluorometric_cdom <= 0) | (ds.fluorometric_cdom > 375)  # fluorometric CDOM measurement range
    cdom_flag[m] = 4
    m = (ds.estimated_chlorophyll <= 0) | (ds.estimated_chlorophyll > 30)  # estimated chlorophyll measurement range
    chl_flag[m] = 4

    return beta_flag, cdom_flag, chl_flag


def flort_datalogger(ds, burst=False):
    """
    Takes FLORT data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial FLORT data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used (from the CSPP data, not part of the mooring data)
    #   measurement_wavelength_* == metadata, move into variable attributes.
    drop_vars = ['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                 'measurement_wavelength_cdom', 'measurement_wavelength_chl']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # check for data from a co-located CTD, if not present add it and reset the fill value for the optical
    # backscatter derived values
    if 'sea_water_temperature' not in ds.variables:
        ds['sea_water_temperature'] = ('time', ds['deployment'].data * FILL_FLOAT)
        ds['sea_water_practical_salinity'] = ('time', ds['deployment'].data * FILL_FLOAT)
        ds['optical_backscatter'] = ds['optical_backscatter'] * FILL_FLOAT
        ds['seawater_scattering_coefficient'] = ds['seawater_scattering_coefficient'] * FILL_FLOAT

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'raw_signal_chl': 'raw_chlorophyll',
        'fluorometric_chlorophyll_a': 'estimated_chlorophyll',
        'fluorometric_chlorophyll_a_qc_executed': 'estimated_chlorophyll_qc_executed',
        'fluorometric_chlorophyll_a_qc_results': 'estimated_chlorophyll_qc_results',
        'raw_signal_cdom': 'raw_fluorometric_cdom',
        'raw_signal_beta': 'raw_backscatter',
        'total_volume_scattering_coefficient': 'beta_700',
        'total_volume_scattering_coefficient_qc_executed': 'beta_700_qc_executed',
        'total_volume_scattering_coefficient_qc_results': 'beta_700_qc_results',
        'optical_backscatter': 'bback',
        'optical_backscatter_qc_executed': 'bback_qc_executed',
        'optical_backscatter_qc_results': 'bback_qc_results',
        'seawater_scattering_coefficient': 'sea_water_scattering_coefficient',
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

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    # check if the older QC flags are present, if not add summary flags with a default value of 1
    flags = ['beta_700_qc_summary_flag', 'fluorometric_cdom_qc_summary_flag', 'estimated_chlorophyll_qc_summary_flag']
    for flag in flags:
        if flag not in ds.variables:
            ds[flag] = ds['time'].astype('int32') * 0 + 1  # default flag values, no errors

    # create QC flags for the data and add them to the QC summary flags
    beta_flag, cdom_flag, chl_flag = quality_checks(ds)
    ds['beta_700_qc_summary_flag'] = ('time', (np.array([ds.beta_700_qc_summary_flag,
                                                         beta_flag])).max(axis=0, initial=1))
    ds['fluorometric_cdom_qc_summary_flag'] = ('time', (np.array([ds.fluorometric_cdom_qc_summary_flag,
                                                                 cdom_flag])).max(axis=0, initial=1))
    ds['estimated_chlorophyll_qc_summary_flag'] = ('time', (np.array([ds.estimated_chlorophyll_qc_summary_flag,
                                                                      chl_flag])).max(axis=0, initial=1))

    if burst:
        # resample the data to the defined time interval
        ds['time'] = ds['time'] + np.timedelta64(450, 's')
        burst = ds.resample(time='900s', skipna=True).median(dim='time', keep_attrs=True)

        # for each of the three FLORT measurements, calculate stats (min, max, and the standard deviation)
        # for each of the bursts
        cdom = ds['fluorometric_cdom'].resample(time='900s', skipna=True)
        cdom = np.array([cdom.min('time').values, cdom.max('time').values, cdom.std('time').values])

        chl = ds['estimated_chlorophyll'].resample(time='900s', skipna=True)
        chl = np.array([chl.min('time').values, chl.max('time').values, chl.std('time').values])

        beta = ds['beta_700'].resample(time='900s', skipna=True)
        beta = np.array([beta.min('time').values, beta.max('time').values, beta.std('time').values])

        # create a data set with the burst statistics for the variables
        stats = xr.Dataset({
            'fluorometric_cdom_burst_stats': (['time', 'stats'], cdom.T),
            'estimated_chlorophyll_burst_stats': (['time', 'stats'], chl.T),
            'beta_700_burst_stats': (['time', 'stats'], beta.T)
        }, coords={'time': burst['time'], 'stats': np.arange(0, 3).astype('int32')})

        # add the stats into the burst averaged data set, and then remove the missing rows
        burst = burst.merge(stats)
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # save the newly average data
        ds = burst

    return ds


def flort_instrument(ds):
    """
    Takes FLORT data recorded by the Sea-Bird Electronics SBE16Plus used in the
    CGSN/EA moorings and cleans up the data set to make it more user-friendly.
    Primary task is renaming parameters and dropping some that are of limited
    use. Additionally, re-organize some of the variables to permit better
    assessments of the data.

    :param ds: initial FLORT data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used (from the CSPP data, not part of the mooring data)
    #   measurement_wavelength_* == metadata, move into variable attributes.
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl'])

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'raw_signal_chl': 'raw_chlorophyll',
        'fluorometric_chlorophyll_a': 'estimated_chlorophyll',
        'fluorometric_chlorophyll_a_qc_executed': 'estimated_chlorophyll_qc_executed',
        'fluorometric_chlorophyll_a_qc_results': 'estimated_chlorophyll_qc_results',
        'raw_signal_cdom': 'raw_fluorometric_cdom',
        'raw_signal_beta': 'raw_backscatter',
        'total_volume_scattering_coefficient': 'beta_700',
        'total_volume_scattering_coefficient_qc_executed': 'beta_700_qc_executed',
        'total_volume_scattering_coefficient_qc_results': 'beta_700_qc_results',
        'optical_backscatter': 'bback',
        'optical_backscatter_qc_executed': 'bback_qc_executed',
        'optical_backscatter_qc_results': 'bback_qc_results',
        'seawater_scattering_coefficient': 'sea_water_scattering_coefficient'
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

    # check if the raw data for all three channels is 0, if so the FLORT wasn't talking to the CTD and these are
    # all just fill values that can be removed.
    ds = ds.where(ds['raw_backscatter'] + ds['raw_fluorometric_cdom'] + ds['raw_chlorophyll'] > 0, drop=True)
    if len(ds.time) == 0:
        # this was one of those deployments where the FLORT was never able to communicate with the CTD.
        warnings.warn('Communication failure between the FLORT and the CTDBP. No data was recorded.')
        return None

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


def flort_cspp(ds):
    """
    Takes FLORT data recorded by the CSPP loggers used by the Endurance Array
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial FLORT data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used, nor well-defined by the CSPP software
    #   measurement_wavelength_* == metadata, move into variable attributes.
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl'])

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'pressure': 'sea_water_pressure',
        'pressure_qc_executed': 'sea_water_pressure_qc_executed',
        'pressure_qc_results': 'sea_water_pressure_qc_results',
        'raw_signal_chl': 'raw_chlorophyll',
        'fluorometric_chlorophyll_a': 'estimated_chlorophyll',
        'fluorometric_chlorophyll_a_qc_executed': 'estimated_chlorophyll_qc_executed',
        'fluorometric_chlorophyll_a_qc_results': 'estimated_chlorophyll_qc_results',
        'raw_signal_cdom': 'raw_fluorometric_cdom',
        'raw_signal_beta': 'raw_backscatter',
        'total_volume_scattering_coefficient': 'beta_700',
        'total_volume_scattering_coefficient_qc_executed': 'beta_700_qc_executed',
        'total_volume_scattering_coefficient_qc_results': 'beta_700_qc_results',
        'optical_backscatter': 'bback',
        'optical_backscatter_qc_executed': 'bback_qc_executed',
        'optical_backscatter_qc_results': 'bback_qc_results',
        'seawater_scattering_coefficient': 'sea_water_scattering_coefficient'
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

    # create a profile variable to uniquely identify profiles within the dataset
    print('Creating and adding a profile variable to the data set ...')
    ds = create_profile_id(ds)

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


def flort_wfp(ds, grid=False):
    """
    Takes FLORT data recorded by the Wire-Following Profilers (used by CGSN/EA
    as part of the coastal and global arrays) and cleans up the data set to
    make it more user-friendly.  Primary task is renaming parameters and
    dropping some that are of limited use. Additionally, re-organize some of
    the variables to permit better assessments of the data.

    :param ds: initial FLORT data set downloaded from OOI via the M2M system
    :param grid: boolean flag for whether the data should be gridded
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used (from the CSPP data, not part of the WFP data)
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   raw_internal_temp == not available, NaN filled
    ds = ds.reset_coords()
    drop_vars = ['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                 'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'raw_internal_temp']
    for v in drop_vars:
        if v in ds.variables:
            ds = ds.drop_vars(v)

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'int_ctd_pressure': 'sea_water_pressure',
        'raw_signal_chl': 'raw_chlorophyll',
        'fluorometric_chlorophyll_a': 'estimated_chlorophyll',
        'fluorometric_chlorophyll_a_qc_executed': 'estimated_chlorophyll_qc_executed',
        'fluorometric_chlorophyll_a_qc_results': 'estimated_chlorophyll_qc_results',
        'raw_signal_cdom': 'raw_fluorometric_cdom',
        'raw_signal_beta': 'raw_backscatter',
        'total_volume_scattering_coefficient': 'beta_700',
        'total_volume_scattering_coefficient_qc_executed': 'beta_700_qc_executed',
        'total_volume_scattering_coefficient_qc_results': 'beta_700_qc_results',
        'optical_backscatter': 'bback',
        'optical_backscatter_qc_executed': 'bback_qc_executed',
        'optical_backscatter_qc_results': 'bback_qc_results',
        'seawater_scattering_coefficient': 'sea_water_scattering_coefficient'
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
                    da.attrs['_FillValue'] = FILL_FLOAT
                    da = da.astype(float)

                # add to the list
                gridded.append(da)

        # recombine the gridded data arrays into a single dataset
        gridded = xr.merge(gridded)
        gridded.attrs = ds.attrs
        ds = gridded
    else:
        # create a profile variable to uniquely identify profiles within the dataset
        print('Creating and adding a profile variable to the data set ...')
        ds = create_profile_id(ds)

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
        flort = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*FLORT.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not flort:
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
        flort = m2m_collect(r, '.*FLORT.*\\.nc$')

        # check to see if we downloaded any data
        if not flort:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize the data
    if node == 'SP001':
        # this FLORT is part of a CSPP
        flort = flort_cspp(flort)
    elif node == 'WFP01':
        # this FLORT is part of a Wire-Following Profiler
        flort = flort_wfp(flort)
    elif node == 'SBD17':
        # this FLORT is connected to the CTDBP on an EA Inshore Surface Mooring
        flort = flort_instrument(flort)
        if not flort:
            # there was no data after removing all the 0's
            sys.exit()
    else:
        # this FLORT is stand-alone on one of the moorings
        flort = flort_datalogger(flort, burst)

    vocab = get_vocabulary(site, node, sensor)[0]
    flort = update_dataset(flort, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    flort.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
