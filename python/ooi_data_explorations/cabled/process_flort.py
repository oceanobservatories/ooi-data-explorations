#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import sys
import time
import xarray as xr

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from tqdm import tqdm

from ooi_data_explorations.common import inputs, load_gc_thredds, m2m_collect, m2m_request, get_vocabulary, \
    update_dataset, ENCODINGS, FILL_FLOAT, N_CORES
from ooi_data_explorations.qartod.qc_processing import parse_qc
from ooi_data_explorations.profilers import updown, create_profile_id, bin_profiles
from ooi_data_explorations.uncabled.process_flort import ATTRS, quality_checks


def flort_profiler(ds):
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
    #   internal_timestamp == not used, nor accurate
    #   time_string == not used, nor accurate (used to construct internal_timestamp)
    #   date_string == not used, nor accurate (used to construct internal_timestamp)
    #   measurement_wavelength_* == metadata, move into variable attributes.
    #   signal_*_scale_factor == filled with 0's, not used
    #   signal_*_offset == filled with 0's, not used
    #   pressure_depth == variable assigned if this was a FLORT on a CSPP, not with moorings
    drop_vars = ['internal_timestamp', 'time_string', 'date_string', 'measurement_wavelength_beta', 
                 'measurement_wavelength_cdom', 'measurement_wavelength_chl', 'signal_1_scale_factor', 
                 'signal_1_offset', 'signal_2_scale_factor', 'signal_2_offset', 'signal_3_scale_factor', 
                 'signal_3_offset']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # check for data from a co-located CTD, if not present add it and reset the fill value for the optical
    # backscatter derived values
    if "sea_water_temperature" not in ds.variables:
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

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    # create QC flags for the data and add them to the OOI QC summary flags
    beta_flag, cdom_flag, chl_flag = quality_checks(ds)
    ds['beta_700_qc_summary_flag'] = ('time', (np.array([ds.beta_700_qc_summary_flag,
                                                         beta_flag])).max(axis=0, initial=1))
    ds['fluorometric_cdom_qc_summary_flag'] = ('time', (np.array([ds.fluorometric_cdom_qc_summary_flag,
                                                                 cdom_flag])).max(axis=0, initial=1))
    ds['estimated_chlorophyll_qc_summary_flag'] = ('time', (np.array([ds.estimated_chlorophyll_qc_summary_flag,
                                                                      chl_flag])).max(axis=0, initial=1))

    # return the re-processed data
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

    flort = flort_profiler(flort)

    vocab = get_vocabulary(site, node, sensor)[0]
    flort = update_dataset(flort, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    flort.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
