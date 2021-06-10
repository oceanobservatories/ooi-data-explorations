#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the PHSEN data from the uncabled, Coastal Endurance Surface
    Moorings and processes the data to generate QARTOD Gross Range and
    Climatology test limits
"""
import datetime
import numpy as np
import os
import pandas as pd
import re
import requests
import sys
import xarray as xr

from ooi_data_explorations.common import gc_collect, add_annotation_qc_flags
from ooi_data_explorations.uncabled.process_phsen import phsen_datalogger, phsen_instrument, quality_checks
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology


def load_gc_thredds(site, node, sensor, method, stream):
    """

    :param site:
    :param node:
    :param sensor:
    :param method:
    :param stream:
    :return:
    """
    # download the data from the Gold Copy THREDDS server
    dataset_id = '-'.join([site, node, sensor, method, stream]) + '/catalog.html'
    tag = '.*PHSEN.*\\.nc$'
    data = gc_collect(dataset_id, tag)

    return data


def combine_datasets():
    pass


def generate_qartod():
    pass


def main():
    pass


if __name__ == '__main__':
    main()
