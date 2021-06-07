#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Calculates QARTOD test ranges and creates the resulting tables used by
    the OOI QC lookup functions implementing the QARTOD testing.
"""
import argparse
import glob
import numpy as np
import os
import sys
import xarray as xr


def calc_gross_range(ds, params):
    """
    Simple function to calculate the user gross range values based on the
    long-term average of the parameter(s) of interest plus/minus 3 standard
    deviations (the statistical rule of three). The user range is used in the
    QARTOD Gross Range test to set a suspect flag for data points that lie
    outside of the defined range.

    :param ds: xarray dataset with the parameter(s) of interest
    :param params: list of the parameter names to calculate the ranges for
    :return user_range: list of gross range values
    """
    user_range = []
    for param in params:
        mu = ds[param].mean().values
        sd = ds[param].std().values * 3
        user_range.append([mu - sd, mu + sd])

    return user_range


def format_gross_range(params, sensor_range, user_range, **kwargs):
    """
    Creates a JSON formatted object that can later be saved to a CSV formatted
    file for use in the Gross Range lookup tables.

    :param params: list of the parameter name(s) corresponding to the
        calculated user ranges
    :param sensor_range: default sensor, or fail range, usually referenced
        from the vendor documentation
    :param user_range: user range, or sensor range, calculated from the data
    :param kwargs: Takes the following optional keyword arguments:

        site: Site designator, extracted from the first part of the reference
            designator
        node: Node designator, extracted from the second part of the reference
            designator
        sensor: Sensor designator, extracted from the third and fourth part of
            the reference designator
        stream: Stream name that contains the data of interest
        source: Notes or comments about how the Gross Range values were
            obtained

    :return gross_range: JSON formatted string with the sensor and user gross
        range values added in the formatting expected by the QC lookup
    """
    pass

