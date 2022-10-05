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
    update_dataset, ENCODINGS
from ooi_data_explorations.qartod.qc_processing import parse_qc

# load configuration settings
FILL_INT = -9999999
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
        'ancillary_variables': ('beta_700 seawater_temperature practical_salinity seawater_scattering_coefficient '
                                'bback_qc_executed bback_qc_results'),
        '_FillValue': np.nan
    },
    'seawater_scattering_coefficient': {
        'long_name': 'Seawater Optical Backscatter at 700 nm',
        'units': 'm-1',
        'comment': ('Theoretical estimation of the optical backscatter for pure seawater at 700 nm adjusted for the'
                    'effects of temperature and salinity. This value is added to the particulate optical backscatter '
                    'measurement to create the total optical backscatter measurement contained in this data set.'),
        'ancillary_variables': 'seawater_temperature practical_salinity',
        '_FillValue': np.nan
    },
    'practical_salinity': {
        'long_name': 'Practical Salinity',
        'standard_name': 'sea_water_practical_salinity',
        'units': '1',
        'comment': ('Salinity is generally defined as the concentration of dissolved salt in a parcel of sea water. ' 
                    'Practical Salinity is a more specific unitless quantity calculated from the conductivity of ' 
                    'sea water and adjusted for temperature and pressure. It is approximately equivalent to Absolute ' 
                    'Salinity (the mass fraction of dissolved salt in sea water), but they are not interchangeable. '
                    'Measurements are from a co-located CTD.'),
        'data_product_identifier': 'PRACSAL_L2',
        '_FillValue': np.nan
    },
    'seawater_temperature': {
        'long_name': 'Seawater Temperature',
        'standard_name': 'sea_water_temperature',
        'units': 'degree_Celsius',
        'comment': ('Sea water temperature is the in situ temperature of the sea water. Measurements are from a '
                    'co-located CTD'),
        'data_product_identifier': 'TEMPWAT_L1',
        '_FillValue': np.nan
    },
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
    max_counts = 4115   # counts should be greater than 0 and less than 4120 +/- 5
    beta_flag = ds['time'].astype('int32') * 0 + 1   # default flag values, no errors
    chl_flag = ds['time'].astype('int32') * 0 + 1   # default flag values, no errors

    # test the min/max values of the raw measurements
    if "raw_backscatter" in ds.variables:
        m = (ds.raw_backscatter <= 0) | (ds.raw_backscatter > max_counts)
        beta_flag[m] = 4    # raw volume scattering coefficient values off scale
    if "raw_chlorophyll" in ds.variables:
        m = (ds.raw_chlorophyll <= 0) | (ds.raw_chlorophyll > max_counts)
        chl_flag[m] = 4     # raw chlorophyll values off scale

    # test the min/max values of the derived measurements (values from the vendor documentation)
    m = (ds.estimated_chlorophyll <= 0) | (ds.estimated_chlorophyll > 30)  # estimated chlorophyll measurement range
    chl_flag[m] = 4

    return beta_flag, chl_flag


def flord_datalogger(ds, burst=False):
    """
    Takes flord data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial flort data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   pressure_depth == variable assigned if this was a FLORT on a CSPP, not with moorings
    drop_vars = ['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                 'measurement_wavelength_cdom', 'measurement_wavelength_chl']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # check for data from a co-located CTD, if not present add it and reset the fill value for the optical
    # backscatter derived values
    if 'temp' not in ds.variables and "seawater_temperature" not in ds.variables:
        ds['temp'] = ('time', ds['deployment'].data * np.nan)
        ds['practical_salinity'] = ('time', ds['deployment'].data * np.nan)

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'temp': 'seawater_temperature',
        'raw_signal_chl_volts': 'raw_chlorophyll_volts',
        'fluorometric_chlorophyll_a': 'estimated_chlorophyll',
        'fluorometric_chlorophyll_a_qc_executed': 'estimated_chlorophyll_qc_executed',
        'fluorometric_chlorophyll_a_qc_results': 'estimated_chlorophyll_qc_results',
        'raw_signal_beta_volts': 'raw_backscatter_volts',
        'total_volume_scattering_coefficient': 'beta_700',
        'total_volume_scattering_coefficient_qc_executed': 'beta_700_qc_executed',
        'total_volume_scattering_coefficient_qc_results': 'beta_700_qc_results',
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

    # create QC flags for the data and add them to the OOI QC summary flags
    beta_flag, chl_flag = quality_checks(ds)
    ds['beta_700_qc_summary_flag'] = ('time', (np.array([ds.beta_700_qc_summary_flag,
                                                         beta_flag])).max(axis=0, initial=1))
    ds['estimated_chlorophyll_qc_summary_flag'] = ('time', (np.array([ds.estimated_chlorophyll_qc_summary_flag,
                                                                      chl_flag])).max(axis=0, initial=1))

    return ds


def flord_instrument(ds):
    """
    Takes flort data recorded by the Sea-Bird Electronics SBE16Plus used in the
    CGSN/EA moorings and cleans up the data set to make it more user-friendly.
    Primary task is renaming parameters and dropping some that are of limited
    use. Additionally, re-organize some of the variables to permit better
    assessments of the data.

    :param ds: initial flort data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   pressure_depth == variable assigned if this was a FLORT on a CSPP, not with moorings
    ds = ds.reset_coords()
    drop_vars = ['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                 'measurement_wavelength_cdom', 'measurement_wavelength_chl']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)
    
    if 'temp' not in ds.variables and "seawater_temperature" not in ds.variables:
        ds['temp'] = ('time', ds['deployment'].data * np.nan)
        ds['practical_salinity'] = ('time', ds['deployment'].data * np.nan)

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'temp': 'seawater_temperature',
        'raw_signal_chl': 'raw_chlorophyll',
        'raw_signal_chl_volts': 'raw_chlorophyll_volts',
        'fluorometric_chlorophyll_a': 'estimated_chlorophyll',
        'fluorometric_chlorophyll_a_qc_executed': 'estimated_chlorophyll_qc_executed',
        'fluorometric_chlorophyll_a_qc_results': 'estimated_chlorophyll_qc_results',
        'raw_signal_beta': 'raw_backscatter',
        'raw_signal_beta_volts': 'raw_backscatter_volts',
        'total_volume_scattering_coefficient': 'beta_700',
        'total_volume_scattering_coefficient_qc_executed': 'beta_700_qc_executed',
        'total_volume_scattering_coefficient_qc_results': 'beta_700_qc_results',
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
    ds = ds.where(ds['raw_backscatter'] + ds['raw_chlorophyll'] > 0, drop=True)
    if len(ds.time) == 0:
        # this was one of those deployments where the FLORT was never able to communicate with the CTD.
        warnings.warn('Communication failure between the FLORD and the CTDBP. No data was recorded.')
        return None

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    # create qc flags for the data and add them to the OOI qc flags
    beta_flag, chl_flag = quality_checks(ds)
    ds['beta_700_qc_summary_flag'] = ('time', (np.array([ds.beta_700_qc_summary_flag,
                                                         beta_flag])).max(axis=0, initial=1))
    ds['estimated_chlorophyll_qc_summary_flag'] = ('time', (np.array([ds.estimated_chlorophyll_qc_summary_flag,
                                                                      chl_flag])).max(axis=0, initial=1))

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
