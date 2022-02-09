#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the FLORT data from the uncabled, Coastal Endurance Surface
    Moorings and Profilers and process the data to generate QARTOD Gross Range
    and Climatology test limits
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz
import xarray as xr

from ooi_data_explorations.common import get_annotations, load_gc_thredds, add_annotation_qc_flags
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_flort import flort_datalogger, flort_instrument, flort_cspp, flort_wfp
from ooi_data_explorations.qartod.qc_processing import identify_blocks, create_annotations, process_gross_range, \
    process_climatology, inputs, ANNO_HEADER, CLM_HEADER, GR_HEADER


def combine_delivery_methods(site, node, sensor):
    """
    Takes the downloaded data from the different data delivery methods for the
    three-channel fluorometer (FLORT), and combines them, where appropriate,
    into a single, merged xarray data sets.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return merged: the merged and resampled (if appropriate) FLORT dataset
    """
    # set the stream and tag constants
    tag = '.*FLORT.*\\.nc$'
    stream = 'flort_sample'

    if node in ['SP001', 'WFP01']:
        # this FLORT is part of a CSPP or WFP and includes telemetered and recovered data
        if node == 'SP001':
            telem = None  # don't use the telemetered CSPP data
            print('##### Downloading the recovered_cspp FLORT data for %s #####' % site)
            rhost = load_gc_thredds(site, node, sensor, 'recovered_cspp', stream, tag)
            deployments = []
            print('# -- Group the data by deployment and process the data')
            grps = list(rhost.groupby('deployment'))
            for grp in grps:
                print('# -- Processing recovered_host deployment %s' % grp[0])
                deployments.append(flort_cspp(grp[1]))
            rhost = xr.concat(deployments, 'time')
        else:
            print('##### Downloading the telemetered FLORT data for %s #####' % site)
            telem = load_gc_thredds(site, node, sensor, 'telemetered', stream, tag)
            deployments = []
            print('# -- Group the data by deployment and process the data')
            grps = list(telem.groupby('deployment'))
            for grp in grps:
                print('# -- Processing telemetered deployment %s' % grp[0])
                deployments.append(flort_wfp(grp[1]))
            telem = xr.concat(deployments, 'time')

            print('##### Downloading the recovered_wfp FLORT data for %s #####' % site)
            rhost = load_gc_thredds(site, node, sensor, 'recovered_wfp', stream, tag)
            deployments = []
            print('# -- Group the data by deployment and process the data')
            grps = list(rhost.groupby('deployment'))
            for grp in grps:
                print('# -- Processing recovered_host deployment %s' % grp[0])
                deployments.append(flort_wfp(grp[1]))
            rhost = xr.concat(deployments, 'time')

        # merge, but do not resample the time records.
        merged = combine_datasets(telem, rhost, None, None)
    elif node == 'SBD17':
        # this FLORT is mounted on the buoy of the Inshore moorings and includes all three types of data
        print('##### Downloading the telemetered FLORT data for %s #####' % site)
        telem = load_gc_thredds(site, node, sensor, 'telemetered', stream, tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(telem.groupby('deployment'))
        for grp in grps:
            print('# -- Processing telemetered deployment %s' % grp[0])
            deployments.append(flort_instrument(grp[1]))
        telem = xr.concat(deployments, 'time')

        print('##### Downloading the recovered_host FLORT data for %s #####' % site)
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host', stream, tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(rhost.groupby('deployment'))
        for grp in grps:
            print('# -- Processing recovered_host deployment %s' % grp[0])
            deployments.append(flort_instrument(grp[1]))
        rhost = xr.concat(deployments, 'time')

        print('##### Downloading the recovered_inst FLORT data for %s #####' % site)
        rinst = load_gc_thredds(site, node, sensor, 'recovered_inst', stream, tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(rinst.groupby('deployment'))
        for grp in grps:
            print('# -- Processing recovered_inst deployment %s' % grp[0])
            deployments.append(flort_instrument(grp[1]))
        rinst = xr.concat(deployments, 'time')

        # merge and resample to a 2 hour data record
        merged = combine_datasets(telem, rhost, rinst, 120)
    else:
        # this FLORT is standalone on one of the NSIFs and includes the telemetered and recovered_host data
        # data is collected in bursts (3 minutes at 1 Hz). process each data set per-deployment
        print('##### Downloading the telemetered FLORT data for %s #####' % site)
        telem = load_gc_thredds(site, node, sensor, 'telemetered', stream, tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(telem.groupby('deployment'))
        for grp in grps:
            print('# -- Processing telemetered deployment %s' % grp[0])
            deployments.append(flort_datalogger(grp[1]))
        telem = xr.concat(deployments, 'time')

        print('##### Downloading the recovered_host FLORT data for %s #####' % site)
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host', stream, tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(rhost.groupby('deployment'))
        for grp in grps:
            print('# -- Processing recovered_host deployment %s' % grp[0])
            deployments.append(flort_datalogger(grp[1]))
        rhost = xr.concat(deployments, 'time')

        # combine the datasets, leaving them as 15-minute median averaged datasets
        merged = combine_datasets(telem, rhost, None, None)

    return merged


def generate_qartod(site, node, sensor, cut_off):
    """
    Load all FLORT data for a defined reference designator (using the site,
    node and sensor names to construct the reference designator) and
    collected via the different data delivery methods and combine them into a
    single data set from which QARTOD test limits for the gross range and
    climatology tests can be calculated.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :param cut_off: string formatted date to use as cut-off for data to add
        to QARTOD test sets
    :return gr_lookup: CSV formatted strings to save to a csv file for the
        QARTOD gross range lookup tables.
    :return clm_lookup: CSV formatted strings to save to a csv file for the
        QARTOD climatology lookup tables.
    :return clm_table: CSV formatted strings to save to a csv file for the
        QARTOD climatology range table for the seafloor pressure and
        temperature.
    """
    # load the combined data for the different sources of FLORT data
    data = combine_delivery_methods(site, node, sensor)

    # create boolean arrays of the data marked as "fail" by the quality checks and generate initial
    # HITL annotations that can be combined with system annotations to create a cleaned up data set
    # prior to calculating the QARTOD test values
    chl_fail = data.estimated_chlorophyll_qc_summary_flag.where(data.estimated_chlorophyll_qc_summary_flag > 3).notnull()
    blocks = identify_blocks(chl_fail, [18, 72])
    chl_hitl = create_annotations(site, node, sensor, blocks)
    chl_hitl['parameters'] = [[22, 1141] for i in chl_hitl['parameters']]

    cdom_fail = data.fluorometric_cdom_qc_summary_flag.where(data.fluorometric_cdom_qc_summary_flag > 3).notnull()
    blocks = identify_blocks(cdom_fail, [18, 72])
    cdom_hitl = create_annotations(site, node, sensor, blocks)
    cdom_hitl['parameters'] = [[23, 1143] for i in cdom_hitl['parameters']]

    beta_fail = data.beta_700_qc_summary_flag.where(data.beta_700_qc_summary_flag > 3).notnull()
    blocks = identify_blocks(beta_fail, [18, 72], 24)
    beta_hitl = create_annotations(site, node, sensor, blocks)
    beta_hitl['parameters'] = [[24, 25, 1139] for i in beta_hitl['parameters']]

    # combine the different dictionaries into a single HITL annotation dictionary for later use
    hitl = chl_hitl.copy()
    for d in (cdom_hitl, beta_hitl):
        for key, value in d.items():
            hitl[key] = hitl[key] + d[key]

    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

    # append the fail annotations to the existing annotations
    annotations = annotations.append(pd.DataFrame(hitl), ignore_index=True, sort=False)

    # create an annotation-based quality flag
    data = add_annotation_qc_flags(data, annotations)

    # clean-up the data, NaN-ing values that were marked as fail in the QC checks and/or identified as a block
    # of failed data, and then removing all records where the rollup annotation (every parameter fails) was
    # set to fail.
    data['estimated_chlorophyll'][chl_fail] = np.nan
    if 'fluorometric_chl_a_annotations_qc_results' in data.variables:
        m = data.fluorometric_chl_a_annotations_qc_results == 4
        data['estimated_chlorophyll'][m] = np.nan

    data['fluorometric_cdom'][cdom_fail] = np.nan
    if 'fluorometric_cdom_annotations_qc_results' in data.variables:
        m = data.fluorometric_cdom_annotations_qc_results == 4
        data['fluorometric_cdom'][m] = np.nan

    data['beta_700'][beta_fail] = np.nan
    if 'total_volume_scattering_coefficient_annotations_qc_results' in data.variables:
        m = data.total_volume_scattering_coefficient_annotations_qc_results == 4
        data['beta_700'][m] = np.nan
        data['bback'][m] = np.nan

    if 'rollup_annotations_qc_results' in data.variables:
        data = data.where(data.rollup_annotations_qc_results < 4)

    # if a cut_off date was used, limit data to all data collected up to the cut_off date.
    # otherwise, set the limit to the range of the downloaded data.
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

    data = data.sel(time=slice('2014-01-01T00:00:00', end_date))

    # set the parameters and the gross range limits
    parameters = ['beta_700', 'bback', 'estimated_chlorophyll', 'fluorometric_cdom']
    limits = [[0, 5], [0, 5], [0, 50], [0, 375]]

    # create the initial gross range entry
    gr_lookup = process_gross_range(data, parameters, limits, site=site, node=node, sensor=sensor)

    # add the stream name and the source comment
    gr_lookup['stream'] = 'flort_sample'
    gr_lookup['source'] = ('Sensor min/max based on the vendor sensor specifications. The user min/max represents '
                           'the range of values that approximately cover 99.7% of all data collected '
                           'through {}.'.format(src_date))

    # create and format the climatology lookups and tables for the data
    clm_lookup, clm_table = process_climatology(data, parameters, limits, site=site, node=node, sensor=sensor)

    # add the stream name
    clm_lookup['stream'] = 'flort_sample'

    return annotations, gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the FLORT data from the Gold Copy THREDDS server and create the
    QARTOD gross range and climatology test lookup tables.
    """
    # setup the input arguments
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    cut_off = args.cut_off

    # create the QARTOD gross range and climatology lookup values and tables
    annotations, gr_lookup, clm_lookup, clm_table = generate_qartod(site, node, sensor, cut_off)

    # save the downloaded annotations and qartod lookups and tables
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/flort')
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
    parameters = ['beta_700', 'bback', 'estimated_chlorophyll', 'fluorometric_cdom']
    for i in range(len(parameters)):
        tbl = '-'.join([site, node, sensor, parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
