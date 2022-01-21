#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import xarray as xr

from ooi_data_explorations.common import inputs, load_gc_thredds, m2m_collect, m2m_request, get_vocabulary, \
    update_dataset, ENCODINGS
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
        'comment': 'Raw chlorophyll fluorescence (470 nm excitation/ 695 nm emission) measurements.',
        'data_product_identifier': 'CHLAFLO_L0'
    },
    'raw_cdom': {
        'long_name': 'Raw CDOM Fluorescence',
        'units': 'counts',
        'comment': 'Raw CDOM fluorescence (370 nm excitation/ 460 nm emission) measurements.',
        'data_product_identifier': 'CDOMFLO_L0'
    },
    'estimated_chlorophyll': {
        'long_name': 'Estimated Chlorophyll Concentration',
        'standard_name': 'mass_concentration_of_chlorophyll_in_sea_water',
        'units': 'ug L-1',
        'comment': ('Estimated chlorophyll concentration based upon a calibration curve derived from a fluorescent ' +
                    'proxy approximately equal to 25 ug/l of a Thalassiosira weissflogii phytoplankton culture. ' +
                    'This measurement is considered to be an estimate only of the true chlorophyll concentration.'),
        'data_product_identifier': 'CHLAFLO_L1',
        'ancillary_variables': 'raw_chlorophyll estimated_chlorophyll_qc_executed estimated_chlorophyll_qc_results'
    },
    'fluorometric_cdom': {
        'long_name': 'Fluorometric CDOM Concentration',
        'standard_name': ('concentration_of_colored_dissolved_organic_matter_in_sea_water_expressed_as_equivalent' +
                          '_mass_fraction_of_quinine_sulfate_dihydrate'),
        'units': 'ppb',
        'comment': ('More commonly referred to as Chromophoric Dissolved Organic Matter (CDOM). CDOM plays an ' +
                    'important role in the carbon cycling and biogeochemistry of coastal waters. It occurs ' +
                    'naturally in aquatic environments primarily as a result of tannins released from decaying ' +
                    'plant and animal matter, and can enter coastal areas in river run-off containing organic ' +
                    'materials leached from soils.'),
        'data_product_identifier': 'CDOMFLO_L1',
        'ancillary_variables': 'raw_cdom fluorometric_cdom_qc_executed fluorometric_cdom_qc_results'
    },
    'beta_700': {
        'long_name': 'Volume Scattering Function at 700 nm',
        'standard_name': 'volume_scattering_function_of_radiative_flux_in_sea_water',
        'units': 'm-1 sr-1',
        'comment': ('Radiative flux is the sum of shortwave and longwave radiative fluxes. Scattering of ' +
                    'radiation is its deflection from its incident path without loss of energy. The volume ' +
                    'scattering function is the intensity (flux per unit solid angle) of scattered radiation per ' +
                    'unit length of scattering medium, normalised by the incident radiation flux.'),
        'data_product_identifier': 'FLUBSCT_L1',
        'ancillary_variables': 'raw_backscatter beta_700_qc_executed beta_700_qc_results'
    },
    'bback': {
        'long_name': 'Total Optical Backscatter at 700 nm',
        'units': 'm-1',
        'comment': ('Total (particulate + water) optical backscatter at 700 nm, derived from the Volume ' +
                    'Scattering Function and corrected for effects of temperature and salinity.'),
        'data_product_identifier': 'FLUBSCT_L2',
        'ancillary_variables': 'beta_700 temperature salinity bback_qc_executed bback_qc_results'
    }
})


def flort_datalogger(ds, burst=True):
    """
    Takes flort data recorded by the data loggers used in the CGSN/EA moorings
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
    #   seawater_scattering_coefficient == not used
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'seawater_scattering_coefficient'])

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

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'temp': 'seawater_temperature',
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

    if burst:   # re-sample the data to a defined time interval using a median average
        # create the burst averaging
        burst = ds
        burst = burst.resample(time='900s', base=3150, loffset='450s', skipna=True).median(dim='time', keep_attrs=True)
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # reset the attributes...which keep_attrs should do...
        burst.attrs = ds.attrs
        for v in burst.variables:
            burst[v].attrs = ds[v].attrs

        # save the newly average data
        ds = burst

    return ds


def flort_instrument(ds):
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
    #   seawater_scattering_coefficient == not used
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'seawater_scattering_coefficient'])

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'temp': 'seawater_temperature',
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
    #   suspect_timestamp = not used
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   seawater_scattering_coefficient == not used
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'seawater_scattering_coefficient'])

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

    return ds


def flort_wfp(ds):
    """
    Takes FLORT data recorded by the Wire-Following Profilers (used by CGSN/EA
    as part of the coastal and global arrays) and cleans up the data set to
    make it more user-friendly.  Primary task is renaming parameters and
    dropping some that are of limited use. Additionally, re-organize some of
    the variables to permit better assessments of the data.

    :param ds: initial FLORT data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = not used
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   seawater_scattering_coefficient == not used
    #   raw_internal_temp == not available, NaN filled
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'seawater_scattering_coefficient',
                  'raw_internal_temp'])

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
    else:
        # this FLORT is on one of the moorings
        if method in ['telemetered', 'recovered_host']:
            flort = flort_datalogger(flort, burst)
        else:
            flort = flort_instrument(flort)

    vocab = get_vocabulary(site, node, sensor)[0]
    flort = update_dataset(flort, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    flort.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
