#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@brief Functions to reprocess and clean up the TURBD datasets downloaded from OOINet
"""

import numpy as np
import pandas as pd
import xarray as xr

from ooi_data_explorations.common import FILL_INT, FILL_FLOAT, ENCODINGS

ATTRS = {
    'raw_backscatter': {
        'long_name': 'Raw Optical Backscatter at 700 nm',
        'units': 'counts',
        'comment': 'Raw optical backscatter measurements at 700 nm.',
        'data_product_identifier': 'FLUBSCT_L0',
        'coordinates': 'time lat lon'
    }
}


def turbd_datalogger(ds, burst=False):
    """
    Takes TURBD data recovred by the data loggers used in CGSN moorings
    and cleeans up the data set to make it more user friendly. Primary 
    task is renaming parameters and dropping some that are of limited use.

    Parameters
    ----------
    ds: xarray.Dataset
        The TURBD dataset with the raw signal and calculated turbidity
    burst: boolean, default=False
        Option to resample the burst sampling of the sensor to median
        average result for the burst duration

    Returns
    -------
    ds: xarray.Dataset
        Cleaned up dataset
    """
    # Drop some redundant variables
    #   internal_timestamp == superceded by time
    #   measurement_wavelength_beta == constant, add as attribute to raw_signal_beta
    drop_vars = ['internal_timestamp', 'measurement_wavelength_beta']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # Rename some variables to have more explanatory names or align with other sensor
    # naming conventions
    rename = {
        'raw_signal_beta': 'raw_backscatter'
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

    # Add in burst resampling
    if burst:
        # resample the data to the defined time interval
        ds['time'] = ds['time'] + np.timedelta64(450, 's')
        burst = ds.resample(time='900s', skipna=True).median(dim='time', keep_attrs=True)

        # for the turbidity, calculate the associated stats for the burst resampling
        turbd = ds['turbidity'].resample(time='900s', skipna=True)
        turbd = np.array([turbd.min('time').values, turbd.max('time').values, turbd.std('time').values])

        # create a data set with the burst statistics for the variables
        stats = xr.Dataset({
            'turbidity_burst_stats': (['time', 'stats'], turbd.T),
        }, coords={'time': burst['time'], 'stats': np.arange(0, 3).astype('int32')})

        # add the stats into the burst averaged data set, and then remove the missing rows
        burst = burst.merge(stats)
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # save the newly average data
        ds = burst

        # Add in attributes for the statistics
        ds['turbidity_burst_stats'].attrs = {
            'long_name': 'Turbidity Statisitics',
            'units': 'ntu',
            'comment': 'Associated statistics for the resampled turbidity measurements.',
            'data_product_identifier': 'FLUBSCT_L0',
            'statistics': 'min max std'
        }

    return ds

