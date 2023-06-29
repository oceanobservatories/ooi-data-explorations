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
from ooi_data_explorations.uncabled.process_flort import ATTRS


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
    :return chl_flag: QARTOD quality flags for the chlorophyll measurements
    """
    max_counts = 4125   # counts should be greater than 0 and less than 4120 +/- 5
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


def flord_instrument(ds):
    """
    Takes FLORD data recorded by the Sea-Bird Electronics SBE16Plus used in the
    CGSN Global moorings and cleans up the data set to make it more user-friendly.
    Primary task is renaming parameters and dropping some that are of limited
    use. Additionally, re-organize some of the variables to permit better
    assessments of the data.

    :param ds: initial flord data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == ctd_time == time, redundant so can remove
    ds = ds.reset_coords()
    drop_vars = ['internal_timestamp', 'ctd_time']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
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

    # check if the raw data for all three channels is 0, if so the FLORD wasn't talking to the CTD and these are
    # all just fill values that can be removed.
    ds = ds.where(ds['raw_backscatter'] + ds['raw_chlorophyll'] > 0, drop=True)
    if len(ds.time) == 0:
        # this was one of those deployments where the FLORD was never able to communicate with the CTD.
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


def flord_wfp(ds, grid=False):
    """
    Takes FLORD data recorded by the Wire-Following Profilers (used by CGSN
    as part of the global arrays) and cleans up the data set to
    make it more user-friendly.  Primary task is renaming parameters and
    dropping some that are of limited use. Additionally, re-organize some of
    the variables to permit better assessments of the data.

    :param ds: initial FLORD data set downloaded from OOI via the M2M system
    :param grid: boolean flag for whether the data should be gridded
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == wfp_timestamp == time, redundant so can remove
    #   raw_internal_temp == not available in WFP data, NaN filled.
    ds = ds.reset_coords()
    drop_vars = ['internal_timestamp', 'wfp_timestamp', 'raw_internal_temp']
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
    beta_flag, chl_flag = quality_checks(ds)
    ds['beta_700_qc_summary_flag'] = ('time', (np.array([ds.beta_700_qc_summary_flag,
                                                         beta_flag])).max(axis=0, initial=1))
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

    # check if we are specifying a deployment or a specific date and time range
    if not deploy or (start and stop):
        return SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')

    # if we are specifying a deployment number, then get the data from the Gold Copy THREDDS server
    if deploy:
        # download the data for the deployment
        flord = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*FLORD.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not flord:
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
        flord = m2m_collect(r, '.*FLORD.*\\.nc$')

        # check to see if we downloaded any data
        if not flord:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize the data
    if node == 'WFP01' or node == 'WFP02':
        # this FLORD is part of a Wire-Following Profiler
        flord = flord_wfp(flord)
    else:
        # this FLORD is connected to a CTDBP
        flord = flord_instrument(flord)
        if not flord:
            # there was no data after removing all the 0's
            sys.exit()

    vocab = get_vocabulary(site, node, sensor)[0]
    flord = update_dataset(flord, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    flord.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
