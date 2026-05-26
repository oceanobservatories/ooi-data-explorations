#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the OPTAA data from the uncabled, Coastal Endurance Surface
    Moorings and Profilers and process the data to generate QARTOD Gross
    Range and Climatology test limits.
"""
# TODO: Methodology for QARTOD test limits is still under development.
# TODO: Processing requires a per-deployment calibration file (cal_file).
#     Determine the standard location and naming convention for cal files
#     before operationalizing this module.
import argparse
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz
import xarray as xr

from ooi_data_explorations.common import get_annotations, get_vocabulary, load_gc_thredds, \
    add_annotation_qc_flags
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_optaa import optaa_datalogger, optaa_cspp
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology, \
    woa_standard_bins, inputs, ANNO_HEADER, CLM_HEADER, GR_HEADER


def combine_delivery_methods(site: str, node: str, sensor: str, cal_file: str) -> xr.Dataset:
    """
    Takes the downloaded data from each of the data delivery methods for
    the OPTAA and combines them into a single, merged xarray data set.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth
        part of the reference designator
    :param cal_file: Path to the AC-S factory calibration file
    :return: merged OPTAA dataset
    """
    tag = '.*OPTAA.*\\.nc$'

    if node == 'SP001':
        # CSPP: recovered_cspp data only
        print('##### Downloading the recovered_cspp OPTAA data for %s #####' % site)
        rcspp = load_gc_thredds(site, node, sensor, 'recovered_cspp',
                                'optaa_dj_cspp_instrument_recovered', tag)
        deployments = []
        grps = list(rcspp.groupby('deployment'))
        for grp in grps:
            print('# -- Processing recovered_cspp deployment %s' % grp[0])
            deployments.append(optaa_cspp(grp[1], cal_file))
        deployments = [i for i in deployments if i]
        rcspp = xr.concat(deployments, 'time')
        merged = combine_datasets(None, rcspp, None, None)
    else:
        # Fixed mooring (NSIF): telemetered and recovered_host
        print('##### Downloading the telemetered OPTAA data for %s #####' % site)
        telem = load_gc_thredds(site, node, sensor, 'telemetered', 'optaa_dj_dcl_instrument', tag)
        deployments = []
        grps = list(telem.groupby('deployment'))
        for grp in grps:
            print('# -- Processing telemetered deployment %s' % grp[0])
            deployments.append(optaa_datalogger(grp[1], cal_file))
        deployments = [i for i in deployments if i]
        telem = xr.concat(deployments, 'time')

        print('##### Downloading the recovered_host OPTAA data for %s #####' % site)
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host',
                                'optaa_dj_dcl_instrument_recovered', tag)
        deployments = []
        grps = list(rhost.groupby('deployment'))
        for grp in grps:
            print('# -- Processing recovered_host deployment %s' % grp[0])
            deployments.append(optaa_datalogger(grp[1], cal_file))
        deployments = [i for i in deployments if i]
        rhost = xr.concat(deployments, 'time')

        merged = combine_datasets(telem, rhost, None, None)

    return merged


def generate_qartod(
    site: str,
    node: str,
    sensor: str,
    cut_off: str,
    cal_file: str
) -> tuple:
    """
    Load all OPTAA data for a defined reference designator and combine
    into a single data set from which QARTOD test limits for the gross
    range and climatology tests can be calculated.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth
        part of the reference designator
    :param cut_off: string formatted date to use as cut-off for data to
        add to QARTOD test sets
    :param cal_file: Path to the AC-S factory calibration file
    :return annotations: DataFrame of system annotations for the sensor
    :return gr_lookup: CSV formatted strings for the QARTOD gross range
        lookup tables
    :return clm_lookup: CSV formatted strings for the QARTOD climatology
        lookup tables
    :return clm_table: CSV formatted strings for the QARTOD climatology
        range tables
    """
    # load the combined data for the different delivery methods
    data = combine_delivery_methods(site, node, sensor, cal_file)

    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(
            annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(
            annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

        # create an annotation-based quality flag
        data = add_annotation_qc_flags(data, annotations)

    if 'rollup_annotations_qc_results' in data.variables:
        data = data.where(data.rollup_annotations_qc_results != 4, drop=True)

    # if a cut_off date was used, limit data to all data collected up to
    # the cut_off date; otherwise set the limit to the data range.
    if cut_off:
        cut = parser.parse(cut_off)
        cut = cut.astimezone(pytz.utc)
        end_date = cut.strftime('%Y-%m-%dT%H:%M:%S')
        src_date = cut.strftime('%Y-%m-%d')
    else:
        cut = parser.parse(data.time_coverage_end)
        cut = cut.astimezone(pytz.utc)
        end_date = cut.strftime('%Y-%m-%dT%H:%M:%S')
        src_date = cut.strftime('%Y-%m-%d')

    _, index = np.unique(data['time'], return_index=True)
    data = data.isel(time=index)
    data = data.sel(time=slice('2014-01-01T00:00:00', end_date))

    # TODO: Define parameters and sensor range limits once the methodology
    #     for OPTAA QARTOD test limits has been established.
    parameters = []
    limits = []

    # set up depth bins for profiling platforms (CSPP)
    depth_bins = np.array([])
    if node == 'SP001':
        vocab = get_vocabulary(site, node, sensor)[0]
        max_depth = vocab['maxdepth']
        depth_bins = woa_standard_bins()
        m = depth_bins[:, 1] <= max_depth
        depth_bins = depth_bins[m, :]

    # TODO: Replace placeholder stream name with the correct stream once
    #     parameters are defined.
    stream = 'optaa_dj_replace_me'

    # create the gross range lookup table
    gr_lookup = process_gross_range(data, parameters, limits, site=site,
                                    node=node, sensor=sensor, stream=stream)
    gr_lookup['source'] = ('User range based on data collected through {}.'.format(src_date))

    # create the climatology lookup and tables
    clm_lookup, clm_table = process_climatology(data, parameters, limits,
                                                depth_bins=depth_bins,
                                                site=site, node=node,
                                                sensor=sensor, stream=stream)

    return annotations, gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the OPTAA data from the Gold Copy THREDDS server and create
    the QARTOD gross range and climatology test lookup tables.
    """
    # extend the standard inputs to include the calibration file path
    args = inputs(argv)
    ap = argparse.ArgumentParser(parents=[args], add_help=False)
    ap.add_argument('-cf', '--cal_file', dest='cal_file', type=str, required=True,
                    help='Path to the AC-S factory calibration file (.dev)')
    args = ap.parse_args(argv)

    site = args.site
    node = args.node
    sensor = args.sensor
    cut_off = args.cut_off
    cal_file = args.cal_file

    # create the QARTOD gross range and climatology lookup values and tables
    annotations, gr_lookup, clm_lookup, clm_table = generate_qartod(
        site, node, sensor, cut_off, cal_file)

    # save the downloaded annotations and qartod lookups and tables
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/optaa')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # save the annotations to a csv file for further processing
    anno_csv = '-'.join([site, node, sensor]) + '.quality_annotations.csv'
    annotations.to_csv(os.path.join(out_path, anno_csv), index=False, columns=ANNO_HEADER)

    # save the gross range values to a csv for further processing
    gr_csv = '-'.join([site, node, sensor]) + '.gross_range.csv'
    gr_lookup.to_csv(os.path.join(out_path, gr_csv), index=False, columns=GR_HEADER)

    # save the climatology values and table to a csv for further processing
    clm_csv = '-'.join([site, node, sensor]) + '.climatology.csv'
    clm_lookup.to_csv(os.path.join(out_path, clm_csv), index=False, columns=CLM_HEADER)
    # TODO: Update parameter list to match the finalized methodology
    parameters = []
    for i in range(len(parameters)):
        tbl = '-'.join([site, node, sensor, parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
