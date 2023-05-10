#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief common functions for processing profiler data
"""
import numpy as np
import xarray as xr


def create_profile_id(ds):
    """
    Use a combination of the deployment number and profile sequence within a
    deployment to create a unique profile identifier for each profile in the
    data set.
    """
    # create an initial profile number for each record in the data set
    ds['profile'] = ds['deployment'].astype('int32') * 0

    # group the data set by deployment
    deployments = ds.groupby('deployment')

    # loop over each deployment
    for deploy in deployments:
        # split the deployment into individual profiles
        profiles = split_profiles(deploy[1])
        profile_id = []
        # add the sequential number for each profile, creating a data array object
        for i, p in enumerate(profiles):
            profile_id.append(p['profile'] + i + 1)

        # combine the profile_id's back into a single data array
        profile_id = xr.concat(profile_id, 'time')

        # update the profile number for each profile in the deployment
        ds['profile'].loc[{'time': profile_id['time']}] = profile_id

    # update the attributes for the profile variable
    ds['profile'].attrs = dict({
        'long_name': 'Profile Number',
        'cf_role': 'profile_id',
        'units': 'count',
        'comment': ('Unique identifier for each profile in a deployment, created by separating the data set into '
                    'individual profiles and then adding a sequential number to each profile. The profile number '
                    'resets to 1 at the start of each deployment. ')
    })
    return ds


def split_profiles(ds):
    """
    Split the data set into individual profiles, where each profile is a
    collection of data from a single deployment and profile sequence. The
    resulting data sets are returned in a list.

    :param ds:
    :return: a list of data sets, one for each profile
    """
    # split the data into profiles, assuming at least 120 seconds between profiles
    dt = ds.where(ds['time'].diff('time') > np.timedelta64(120, 's'), drop=True).get_index('time')

    # process each profile, adding the results to a list of profiles
    profiles = []
    jback = np.timedelta64(30, 's')  # 30 second jump back to avoid collecting data from the following profile
    for i, d in enumerate(dt):
        # pull out the profile
        if i == 0:
            profile = ds.sel(time=slice(ds['time'].values[0], d - jback))
        else:
            profile = ds.sel(time=slice(dt[i - 1], d - jback))

        # add the profile to the list
        profiles.append(profile)

    # grab the last profile and append it to the list
    profile = ds.sel(time=slice(d, ds['time'].values[-1]))
    profiles.append(profile)
    return profiles
