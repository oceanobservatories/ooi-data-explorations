#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the PHSEN data from the uncabled, Coastal Endurance Surface
    Moorings and processes the data to generate QARTOD Gross Range and
    Climatology test limits
"""
import os

from ooi_data_explorations.common import list_deployments, get_deployment_dates, get_vocabulary, m2m_request, \
    m2m_collect, update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_phsen import phsen_datalogger, phsen_instrument, quality_checks
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology


def load_datasets():
    pass


def combine_datasets():
    pass


def generate_qartod():
    pass


def main():
    pass


if __name__ == '__main__':
    main()
