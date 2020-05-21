#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, get_deployment_dates, \
    get_vocabulary, dt64_epoch, update_dataset, ENCODINGS

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
    Takes flort data recorded by the data loggers used in the CGSN/EA moorings and cleans up the data set to make
    it more user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that
    are of no use/value.

    :param ds: initial flort data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = meaningless, should never have been created
    #   provenance == better to access with direct call to OOI M2M api, it doesn't work well in this format
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   pressure_depth == variable assigned if this was a flort on a CSPP, but with moorings
    #   seawater_scattering_coefficient == can safely be ignored, expert users can calculate if they want.
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'pressure_depth',
                  'seawater_scattering_coefficient'])

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
            'units': 'degree_Celsius',
            'instrument': 'CE01ISSM-RID16-03-CTDBPC000',
            'stream': 'ctdbp_cdef_dcl_instrument'
        }

        ds['practical_salinity'] = ('time', ds['deployment'] * np.nan)
        ds['practical_salinity'].attrs = {
            'long_name': 'Practical Salinity',
            'standard_name': 'sea_water_practical_salinity',
            'units': '1',
            'comment': ('Normally this would be seawater salinity data from a co-located CTD. However, data from ' +
                        'that sensor is unavailable. This value has been filled with NaNs to preserve the structure ' +
                        'of the data set.'),
            'data_product_identifier': 'PRACSAL_L2',
            'instrument': 'CE01ISSM-RID16-03-CTDBPC000',
            'stream': 'ctdbp_cdef_dcl_instrument'
        }

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'temp': 'temperature',
        'practical_salinity': 'salinity',
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
    for key, value in rename.items():   # bulk attribute update...
        if value in ATTRS.keys():
            ds[value].attrs = ATTRS[value]
        ds[value].attrs['ooinet_variable_name'] = key
    # ...and for cdom
    ds['fluorometric_cdom'].attrs = ATTRS['fluorometric_cdom']

    # correct incorrect units
    ds['temperature'].attrs['units'] = 'degree_Celsius'

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


def flort_instrument(ds):
    """
    Takes flort data recorded by the Sea-Bird Electronics SBE16Plus used in the CGSN/EA moorings and cleans up the
    data set to make it more user-friendly. Primary task is renaming the alphabet soup parameter names and dropping
    some parameters that are of no use/value.

    :param ds: initial flort data set downloaded from OOI via the M2M system
    :return: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == superseded by time, redundant so can remove
    #   suspect_timestamp = meaningless, should never have been created
    #   provenance == better to access with direct call to OOI M2M api, it doesn't work well in this format
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   pressure_depth == variable assigned if this was a flort on a CSPP, but with moorings
    #   seawater_scattering_coefficient == can safely be ignored, expert users can calculate if they want.
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp', 'suspect_timestamp', 'measurement_wavelength_beta',
                  'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'pressure_depth',
                  'seawater_scattering_coefficient'])

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'ctdbp_seawater_temperature': 'temperature',
        'practical_salinity': 'salinity',
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
    for key, value in rename.items():   # bulk attribute update...
        if value in ATTRS.keys():
            ds[value].attrs = ATTRS[value]
        ds[value].attrs['ooinet_variable_name'] = key
    # ...and for cdom
    ds['fluorometric_cdom'].attrs = ATTRS['fluorometric_cdom']

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
        flort = m2m_collect(r, '.*deployment%04d.*FLORT.*\\.nc$')
    else:
        flort = m2m_collect(r, '.*FLORT.*\\.nc$')

    if not flort:
        exit_text = ('Data unavailable for %s-%s-%s. Check request.' % (site, node, sensor))
        raise SystemExit(exit_text)

    # clean-up and reorganize
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
