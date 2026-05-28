#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the FLORT data from the IOOS GliderDAC for the OOI Coastal
    Endurance gliders and process the data to generate QARTOD Gross Range
    and Climatology test limits.
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz

from ooi_data_explorations.common import get_vocabulary
from ooi_data_explorations.gdac import collect_glider, glider_institution, glider_search_for
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology, \
    woa_standard_bins, inputs, CLM_HEADER, GR_HEADER


def combine_delivery_methods(site: str, node: str, sensor: str):
    """
    Collect all CE glider optical data from the GliderDAC near the
    deployment location of the given reference designator.

    Vocabulary is used to determine the latitude and longitude used as
    the center of the GliderDAC search.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and
        fourth part of the reference designator
    :return: xarray Dataset of all CE glider optical backscatter,
        chlorophyll, and CDOM data near the site
    """
    institution = glider_institution(site)
    search_for = glider_search_for(site)

    data = collect_glider(institution, search_for,
                          latitude=[43.0, 48.0], longitude=[-128.0, -123.9],
                          variables=['bback', 'estimated_chlorophyll', 'fluorometric_cdom'])
    return data


def generate_qartod(site: str, node: str, sensor: str, cut_off: str) -> tuple:
    """
    Load all CE glider optical data from the GliderDAC and combine into
    a single data set from which QARTOD test limits for the gross range
    and climatology tests can be calculated.

    GliderDAC data carries no OOI system annotations, so no annotation
    filtering is applied.

    Both telemetered and recovered_host FLORT glider data use the stream
    name ``flort_m_sample``, so a single set of lookup table entries
    covers both delivery methods.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and
        fourth part of the reference designator
    :param cut_off: string formatted date to use as cut-off for data to
        add to QARTOD test sets
    :return annotations: empty DataFrame (no OOI annotations for GDAC
        data)
    :return gr_lookup: CSV formatted strings for the QARTOD gross range
        lookup tables
    :return clm_lookup: CSV formatted strings for the QARTOD
        climatology lookup tables
    :return clm_table: CSV formatted strings for the QARTOD climatology
        range tables
    """
    data = combine_delivery_methods(site, node, sensor)

    if cut_off:
        cut = parser.parse(cut_off)
        cut = cut.astimezone(pytz.utc)
        end_date = cut.strftime('%Y-%m-%dT%H:%M:%S')
        src_date = cut.strftime('%Y-%m-%d')
    else:
        cut = parser.parse(str(data.time.values.max())[:19])
        cut = cut.astimezone(pytz.utc)
        end_date = cut.strftime('%Y-%m-%dT%H:%M:%S')
        src_date = cut.strftime('%Y-%m-%d')

    _, index = np.unique(data['time'], return_index=True)
    data = data.isel(time=index)
    data = data.sel(time=slice('2014-01-01T00:00:00', end_date))

    # set the parameters and sensor range limits
    parameters = ['bback', 'estimated_chlorophyll', 'fluorometric_cdom']
    limits = [[0, 3], [0, 30], [0, 375]]

    # set up depth bins -- gliders profile the water column
    vocab = get_vocabulary(site, node, sensor)[0]
    max_depth = vocab['maxdepth']
    depth_bins = woa_standard_bins()
    m = depth_bins[:, 1] <= max_depth
    depth_bins = depth_bins[m, :]

    # both delivery methods share the same stream name for glider FLORT
    stream = 'flort_m_sample'

    # create the gross range lookup table
    gr_lookup = process_gross_range(data, parameters, limits, site=site, node=node,
                                    sensor=sensor, stream=stream)
    gr_lookup['source'] = ('User range based on data collected through {}.'.format(src_date))

    # create the climatology lookup and tables
    clm_lookup, clm_table = process_climatology(data, parameters, limits, depth_bins=depth_bins,
                                                site=site, node=node, sensor=sensor, stream=stream)
    clm_lookup['source'] = ('Climatology based on data collected through {}.'.format(src_date))

    return gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the FLORT data from the IOOS GliderDAC and create the
    QARTOD gross range and climatology test lookup tables.
    """
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    cut_off = args.cut_off

    gr_lookup, clm_lookup, clm_table = generate_qartod(site, node, sensor, cut_off)

    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/flort_glider')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    gr_csv = '-'.join([site, node, sensor]) + '.gross_range.csv'
    gr_lookup.to_csv(os.path.join(out_path, gr_csv), index=False, columns=GR_HEADER)

    clm_csv = '-'.join([site, node, sensor]) + '.climatology.csv'
    clm_lookup.to_csv(os.path.join(out_path, clm_csv), index=False, columns=CLM_HEADER)

    parameters = ['bback', 'estimated_chlorophyll', 'fluorometric_cdom']
    for i in range(len(parameters)):
        tbl = '-'.join([site, node, sensor, parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
