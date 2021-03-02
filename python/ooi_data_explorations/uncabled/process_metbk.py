#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, get_deployment_dates, \
    get_vocabulary, dt64_epoch, update_dataset, ENCODINGS
from gsw.conversions import SP_from_C


def metbk_hourly(ds):
    """
    Takes METBK hourly averaged bulk flux estimates from the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly. Primary task is
    renaming parameters and dropping some that are of limited use.
    Additionally, re-organize some of the variables to permit better
    assessments of the data.

    :param ds: initial metbk hourly averaged data set downloaded from OOI via
        the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   met_timeflx == time, redundant,
    #   ### Data products from upstream processing used to calculate hourly flux measurements. Remove from here to
    #   ### keep this data set clean. Will obtain the 1 minute source data from a separate stream.
    #   eastward_velocity
    #   northward_velocity
    #   longwave_irradiance
    #   air_temperature
    #   barometric_pressure
    #   precipitation
    #   sea_surface_temperature
    #   relative_humidity
    #   shortwave_irradiance
    ds = ds.drop(['met_timeflx', 'eastward_velocity', 'northward_velocity', 'longwave_irradiance', 'air_temperature',
                  'barometric_pressure', 'precipitation', 'sea_surface_temperature', 'relative_humidity',
                  'shortwave_irradiance'])

    # reset incorrectly formatted temperature units
    temp_vars = ['met_tempa2m', 'met_tempskn']
    for var in temp_vars:
        ds[var].attrs['units'] = 'degree_Celsius'

    return ds


def metbk_datalogger(ds, burst=False):
    """
    Takes METBK data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use.
    Additionally, re-organize some of the variables to permit better
    assessments of the data.

    :param ds: initial metbk data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   date_time_string == internal_timestamp, redundant so can remove
    #   dcl_controller_timestamp == time, redundant so can remove
    #   internal_timestamp == doesn't exist, always empty so can remove
    #   ### Data products from downstream processing used to calculate hourly flux measurements. Remove from here to
    #   ### keep this data set clean. Will obtain hourly flux data from a different stream.
    #   met_barpres
    #   met_windavg_mag_corr_east
    #   met_windavg_mag_corr_north
    #   met_netsirr
    #   met_salsurf
    #   met_spechum
    #   ct_depth
    #   met_current_direction
    #   met_current_speed
    #   met_relwind_direction
    #   met_relwind_speed
    #   met_heatflx_minute
    #   met_latnflx_minute
    #   met_netlirr_minute
    #   met_sensflx_minute
    ds = ds.drop(['dcl_controller_timestamp', 'internal_timestamp', 'met_barpres',
                  'met_windavg_mag_corr_east', 'met_windavg_mag_corr_north', 'met_netsirr', 'met_salsurf',
                  'met_spechum', 'ct_depth', 'met_current_direction', 'met_current_speed', 'met_relwind_direction',
                  'met_relwind_speed', 'met_heatflx_minute', 'met_latnflx_minute', 'met_netlirr_minute',
                  'met_sensflx_minute', 'met_barpres_qc_executed', 'met_barpres_qc_results',
                  'met_current_direction_qc_executed', 'met_current_direction_qc_results',
                  'met_current_speed_qc_executed', 'met_current_speed_qc_results', 'met_relwind_direction_qc_executed',
                  'met_relwind_direction_qc_results', 'met_relwind_speed_qc_executed', 'met_relwind_speed_qc_results',
                  'met_netsirr_qc_executed', 'met_netsirr_qc_results', 'met_salsurf_qc_executed',
                  'met_salsurf_qc_results', 'met_spechum_qc_executed', 'met_spechum_qc_results'])

    # drop the QC test applied to the L0 values (not supposed to happen)
    ds = ds.drop(['precipitation_qc_executed', 'precipitation_qc_results'])

    # reset incorrectly formatted temperature and relative humidity units
    ds['relative_humidity'].attrs['units'] = 'percent'
    temp_vars = ['air_temperature', 'sea_surface_temperature']
    for var in temp_vars:
        ds[var].attrs['units'] = 'degree_Celsius'

    # calculate the near surface salinity
    ds['sea_surface_salinity'] = ('time', SP_from_C(ds['sea_surface_conductivity'] * 10, ds['sea_surface_temperature'],
                                                    1.0))
    ds['sea_surface_salinity'].attrs = {
        'long_name': 'Sea Surface Practical Salinity',
        'standard_name': 'sea_surface_salinity',
        'units': '1e-3',
        'comment': ('Salinity is generally defined as the concentration of dissolved salt in a parcel of seawater. ' +
                    'Practical Salinity is a more specific unitless quantity calculated from the conductivity of ' +
                    'seawater and adjusted for temperature and pressure. It is approximately equivalent to Absolute ' +
                    'Salinity (the mass fraction of dissolved salt in seawater), but they are not interchangeable.'),
        'data_product_identifier': 'SALSURF_L2',
        'instrument': (ds.attrs['subsite'] + '-SBD11-06-METBKA000'),
        'stream': ds.attrs['stream'],
        'ancillary_variables': 'sea_surface_conductivity sea_surface_temperature'
    }

    if burst:   # re-sample the data to a 15 minute interval using a median average
        burst = ds
        burst = burst.resample(time='900s', base=3150, loffset='450s', keep_attrs=True, skipna=True).median()
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # reset the attributes...which keep_attrs should do...
        burst.attrs = ds.attrs
        for v in burst.variables:
            burst[v].attrs = ds[v].attrs

        # save the newly average data
        ds = burst

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
    if stream in ['metbk_a_dcl_instrument', 'metbk_a_dcl_instrument_recovered']:
        if deploy:
            metbk = m2m_collect(r, ('.*deployment%04d.*METBK.*\\.nc$' % deploy))
        else:
            metbk = m2m_collect(r, '.*METBK.*\\.nc$')

        if not metbk:
            exit_text = ('Data unavailable for %s-%s-%s. Check request.' % (site, node, sensor))
            raise SystemExit(exit_text)

        # clean-up and reorganize
        metbk = metbk_datalogger(metbk)

    else:
        if deploy:
            metbk = m2m_collect(r, ('.*deployment%04d.*METBK.*hourly.*\\.nc$' % deploy))
        else:
            metbk = m2m_collect(r, '.*METBK.*hourly.*\\.nc$')

        if not metbk:
            exit_text = ('Data unavailable for %s-%s-%s. Check request.' % (site, node, sensor))
            raise SystemExit(exit_text)

        # clean-up and reorganize
        metbk = metbk_hourly(metbk)

    vocab = get_vocabulary(site, node, sensor)[0]
    metbk = update_dataset(metbk, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    metbk.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
