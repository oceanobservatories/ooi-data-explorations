#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the DOSTA data from the IOOS GliderDAC for the OOI Coastal
    Endurance gliders and process the data to generate QARTOD Gross Range
    and Climatology test limits.
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz

from ooi_data_explorations.common import list_nodes
from ooi_data_explorations.gdac import collect_glider, glider_institution
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology, \
    woa_standard_bins, inputs, CLM_HEADER, GR_HEADER


def combine_delivery_methods(site: str):
    """
    Collect all CE glider dissolved oxygen data from the GliderDAC near
    the deployment location of the given reference designator.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :return: xarray Dataset of all CE glider dissolved oxygen data
        near the site
    """
    institution = glider_institution(site)
    search_for = None  # all CE data sets, delayed and/or realtime

    data = collect_glider(institution, search_for,
                          latitude=[43.0, 48.0], longitude=[-128.0, -123.9],
                          variables=['dissolved_oxygen'], subsample=33)
    return data


def generate_qartod(site: str, cut_off: str) -> tuple:
    """
    Load all CE glider dissolved oxygen data from the GliderDAC and
    combine into a single data set from which QARTOD test limits for
    the gross range and climatology tests can be calculated.

    GliderDAC data carries no OOI system annotations, so no annotation
    filtering is applied.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param cut_off: string formatted date to use as cut-off for data to
        add to QARTOD test sets
    :return gr_lookup: CSV formatted strings for the QARTOD gross range
        lookup tables
    :return clm_lookup: CSV formatted strings for the QARTOD
        climatology lookup tables
    :return clm_table: CSV formatted strings for the QARTOD climatology
        range tables
    """
    data = combine_delivery_methods(site)

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
    parameters = ['dissolved_oxygen']
    limits = [[0, 500]]

    # OOI stream names for both delivery methods
    streams = ['dosta_abcdjm_glider_instrument', 'dosta_abcdjm_glider_recovered']

    # OOI nodes (glider serial numbers) for the Endurance Array
    nodes = list_nodes(site)
    sensor = '04-DOSTAM000'

    # create the gross range lookup table and replicate for all nodes and streams
    gr_lookup = process_gross_range(data, parameters, limits, site=site, node=nodes[0],
                                    sensor=sensor, stream=streams[0])
    gr_lookup = pd.concat([gr_lookup] * 2 * len(nodes), ignore_index=True)

    idx = 0
    for node in nodes:
        for stream in streams:
            for j in range(len(parameters)):
                gr_lookup.loc[idx + j, 'node'] = node
                gr_lookup.loc[idx + j, 'stream'] = stream
            idx += len(parameters)

    gr_lookup['source'] = ('User range based on data collected through {}.'.format(src_date))

    # set up depth bins -- gliders profile the water column
    depth_bins = woa_standard_bins()
    m = depth_bins[:, 1] <= 1000  # if a glider goes past 1000 m, it isn't coming back
    depth_bins = depth_bins[m, :]

    # create the climatology lookup and tables and replicate for all nodes and streams
    clm_lookup, clm_table = process_climatology(data, parameters, limits, depth_bins=depth_bins,
                                                site=site, node=nodes[0], sensor=sensor,
                                                stream=streams[0])
    clm_lookup = pd.concat([clm_lookup] * 2 * len(nodes), ignore_index=True)

    idx = 0
    for node in nodes:
        for stream in streams:
            for j in range(len(parameters)):
                clm_lookup.loc[idx + j, 'node'] = node
                clm_lookup.loc[idx + j, 'stream'] = stream
            idx += len(parameters)

    clm_lookup['source'] = (
        'Climatology based on depth bins (from {} to {} m), '
        'using data collected through {}.'.format(
            int(depth_bins[0, 0]), int(depth_bins[-1, 1]), src_date
        )
    )

    return gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the DOSTA data from the IOOS GliderDAC and create the
    QARTOD gross range and climatology test lookup tables.
    """
    args = inputs(argv)
    site = args.site
    cut_off = args.cut_off

    gr_lookup, clm_lookup, clm_table = generate_qartod(site, cut_off)

    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/glider')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    gr_csv = site.lower() + '.dosta.gross_range.csv'
    gr_lookup.to_csv(os.path.join(out_path, gr_csv), index=False, columns=GR_HEADER)

    clm_csv = site.lower() + '.dosta.climatology.csv'
    clm_lookup.to_csv(os.path.join(out_path, clm_csv), index=False, columns=CLM_HEADER)

    for i in range(len(parameters := ['dissolved_oxygen'])):
        tbl = '.'.join([site.lower(), 'dosta.climatology', parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
