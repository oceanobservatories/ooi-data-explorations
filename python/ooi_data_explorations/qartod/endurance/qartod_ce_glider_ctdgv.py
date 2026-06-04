#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the CTDGV data from the IOOS GliderDAC for the OOI Coastal
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


def combine_delivery_methods(qc_flags: bool = False):
    """
    Collect 33% of all CE glider CTD data from the GliderDAC for the
    Endurance array.

    :param qc_flags: When True, also download the QARTOD primary summary
        flags for each CTD variable (see :func:`collect_glider`)
    :return: xarray Dataset of all CE glider CTD data
    """
    institution = glider_institution('CE05MOAS')
    search_for = None  # all CE data sets, delayed and/or realtime

    data = collect_glider(institution, search_for,
                          latitude=[43.0, 48.0], longitude=[-128.0, -123.9],
                          variables=['conductivity', 'temperature', 'salinity'],
                          qc_flags=qc_flags, subsample=15)
    return data


def generate_qartod(site: str, cut_off: str) -> tuple:
    """
    Load all CE glider CTD data from the GliderDAC and combine into a
    single data set from which QARTOD test limits for the gross range
    and climatology tests can be calculated.

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
    data = combine_delivery_methods()

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
    parameters = ['conductivity', 'temperature', 'pressure', 'salinity']
    limits = [[0, 9], [-5, 35], [0, 1500], [0, 42]]

    # OOI stream names for both delivery methods
    streams = ['ctdgv_m_glider_instrument', 'ctdgv_m_glider_instrument_recovered']

    # OOI nodes (glider serial numbers) for the Endurance Array
    nodes = list_nodes(site)

    # create the gross range lookup table and replicate for both streams
    sensor = '05-CTDGVM000'
    gr_lookup = process_gross_range(data, parameters, limits, site=site, node=nodes[0],
                                    sensor=sensor, stream=streams[0], stdx=5)
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

    # create the climatology lookup and tables -- temperature and salinity only
    parameters = ['temperature', 'salinity']
    limits = [[-5, 35], [0, 42]]
    clm_lookup, clm_table = process_climatology(data, parameters, limits, depth_bins=depth_bins, site=site,
                                                node=nodes[0], sensor=sensor, stream=streams[0], stdx=5)
    clm_lookup = pd.concat([clm_lookup] * 2 * len(nodes), ignore_index=True)

    idx = 0
    for node in nodes:
        for stream in streams:
            for j in range(2):
                clm_lookup.loc[idx + j, 'node'] = node
                clm_lookup.loc[idx + j, 'stream'] = stream
            idx += 2

    clm_lookup['source'] = (
        'Climatology based on depth bins (from {} to {} m), '
        'using data collected through {}.'.format(
            int(depth_bins[0, 0]), int(depth_bins[-1, 1]), src_date
        )
    )

    return gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the CTDGV data from the IOOS GliderDAC and create the
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

    gr_csv = site.lower() + '.ctdgv.gross_range.csv'
    gr_lookup.to_csv(os.path.join(out_path, gr_csv), index=False, columns=GR_HEADER)

    clm_csv = site.lower() + '.ctdgv.climatology.csv'
    clm_lookup.to_csv(os.path.join(out_path, clm_csv), index=False, columns=CLM_HEADER)

    parameters = ['temperature', 'salinity']
    for i in range(len(parameters)):
        tbl = '.'.join([site.lower(), 'ctdgv.climatology', parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
