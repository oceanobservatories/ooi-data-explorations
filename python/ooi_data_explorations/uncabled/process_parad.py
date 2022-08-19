#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import pandas as pd
import xarray as xr

from scipy.interpolate import griddata

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, load_gc_thredds, \
    get_vocabulary, update_dataset, ENCODINGS
from ooi_data_explorations.qartod.qc_processing import parse_qc

FILL_INT = -9999999
ATTRS = {
    'raw_par_measurement': {
        'long_name': 'Raw PAR Measurement',
        'units': '',  # deliberately left blank, will fill in below
        'comment': 'Photosynthetically Active Radiation (PAR) unprocessed sensor reading.',
        'data_product_identifier': 'OPTPARW_L0'
    },
    'par': {
        'long_name': 'Photosynthetically Active Radiation',
        'standard_name': 'downwelling_photosynthetic_photon_flux_in_sea_water',
        'units': 'umol m-2 s-1',
        'comment': ('Photosynthetically Active Radiation (PAR), or photosynthetic photon flux density, is a measure '
                    'of the density of photons per unit area that are in the spectral range of light (400-700 '
                    'nanometers) used by primary producers photosynthesis.'),
        'data_product_identifier': 'OPTPARW_L1',
        'ancillary_variables': 'raw_par_measurement'
    },
    'seawater_pressure': {
        'long_name': 'Sea Water Pressure',
        'standard_name': 'sea_water_pressure_due_to_sea_water',
        'units': 'dbar',
        'comment': ('Sea Water Pressure refers to the pressure exerted on a sensor in situ by the weight of the ' 
                    'column of seawater above it. It is calculated by subtracting one standard atmosphere from the ' 
                    'absolute pressure at the sensor to remove the weight of the atmosphere on top of the water ' 
                    'column. The pressure at a sensor in situ provides a metric of the depth of that sensor. '
                    'Measurements are from a co-located CTD.'),
        'data_product_identifier': 'PRESWAT_L1'
    }
}


def parad_cspp(ds):
    """
    Takes PARAD data recorded by the CSPP loggers used by the Endurance Array
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some variables to permit better assessments of the data.

    :param ds: initial PARAD data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   suspect_timestamp = not used
    #   internal_timestamp = not accurate, don't use
    #   profiler_timestamp == time, redundant
    #   date_string == used to set internal_timestamp, redundant and not accurate
    #   time_string == used to set internal_timestamp, redundant and not accurate
    drop_list = ['suspect_timestamp', 'internal_timestamp', 'profiler_timestamp', 'date_string', 'time_string']
    for var in ds.variables:
        if var in drop_list:
            ds = ds.drop(var)

    # rename some variables for better clarity
    rename = {
        'par': 'raw_par_measurement',
        'parad_j_par_counts_output': 'par',
        'parad_j_par_counts_output_qc_executed': 'par_qc_executed',
        'parad_j_par_counts_output_qc_results': 'par_qc_results',
        'pressure': 'seawater_pressure',
        'pressure_qc_executed': 'seawater_pressure_qc_executed',
        'pressure_qc_results': 'seawater_pressure_qc_results',
    }
    ds = ds.rename(rename)

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # add original OOINet variable name as an attribute if renamed
    for key, value in rename.items():
        ds[value].attrs['ooinet_variable_name'] = key

    # set the units for the raw PAR measurement
    ds['raw_par_measurement'].attrs['units'] = 'count'

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    return ds


def parad_wfp(ds, grid=False):
    """
    Takes PARAD data recorded by the Wire-Following Profilers (used by CGSN/EA
    as part of the coastal and global arrays) and cleans up the data set to
    make it more user-friendly.  Primary task is renaming parameters and
    dropping some that are of limited use. Additionally, re-organize some of
    the variables to permit better assessments of the data.

    :param ds: initial PARAD data set downloaded from OOI via the M2M system
    :param grid: boolean flag for whether the data should be gridded
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == time, redundant so can remove
    ds = ds.reset_coords()
    ds = ds.drop(['internal_timestamp'])

    # lots of renaming here to get a better defined data set with cleaner attributes
    rename = {
        'int_ctd_pressure': 'seawater_pressure',
        'par_val_v': 'raw_par_measurement',
        'parad_k_par': 'par',
        'parad_k_par_qc_executed': 'par_qc_executed',
        'parad_k_par_qc_results': 'par_qc_results'
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

    # set the units for the raw PAR measurement
    ds['raw_par_measurement'].attrs['units'] = 'mV'

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

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
        # the record to the end.
        depth_range = np.arange(30, 511, 1)
        time_range = np.arange(0, np.ceil(dt.max()) + 1, 8/24)
        gridded_time = base_time.astype('M8[D]') + pd.to_timedelta(time_range, unit='D')

        # grid the data, adding the results to a list of data arrays
        gridded = []
        for v in ds.variables:
            if v not in ['time', 'depth']:
                # grid the data for each variable
                gdata = griddata((dt.values, ds['depth'].values), ds[v].values,
                                 (time_range[None, :], depth_range[:, None]),
                                 method='nearest')

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
    grid = args.burst

    # check if we are specifying a deployment or a specific date and time range
    if not deploy or (start and stop):
        return SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')

    # if we are specifying a deployment number, then get the data from the Gold Copy THREDDS server
    if deploy:
        # download the data for the deployment
        parad = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*PARAD.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not parad:
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
        parad = m2m_collect(r, '.*PARAD.*\\.nc$')

        # check to see if we downloaded any data
        if not parad:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize
    if node == 'SP001':
        parad = parad_cspp(parad)
    elif node == 'WFP01':
        parad = parad_wfp(parad, grid=grid)
    else:
        exit_text = ('Unknown PAR sensor, the node name %s is unrecognized.' % node)
        raise SystemExit(exit_text)

    vocab = get_vocabulary(site, node, sensor)[0]
    parad = update_dataset(parad, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    parad.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
