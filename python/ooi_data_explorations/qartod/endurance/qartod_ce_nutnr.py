#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the NUTNR data from the uncabled Coastal Endurance Surface Moorings
    and process the data to generate QARTOD Gross Range and Climatology test
    limits
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz
import xarray as xr

from ooi_data_explorations.common import get_annotations, get_vocabulary, load_gc_thredds, add_annotation_qc_flags
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_nutnr import suna_datalogger, suna_instrument
from ooi_data_explorations.qartod.qc_processing import identify_blocks, create_annotations, process_gross_range, \
    process_climatology, woa_standard_bins, inputs, ANNO_HEADER, CLM_HEADER, GR_HEADER


def combine_delivery_methods(site, node, sensor):
    """
    Takes the downloaded data from the different data delivery methods for the
    SUNA nitrate sensor (NUTNR), and combines them, where appropriate,
    into a single, merged xarray data sets.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return merged: the merged and resampled (if appropriate) NUTNR dataset
    """
    # set the stream and tag constants
    tag = '.*NUTNR.*\\.nc$'
    stream = ['suna_dcl_recovered', 'suna_dcl_recovered', 'suna_instrument_recovered']

    print('##### Downloading the telemetered NUTNR data for %s #####' % site)
    telem = load_gc_thredds(site, node, sensor, 'telemetered', stream[0], tag)
    deployments = []
    print('# -- Group the data by deployment and process the data')
    grps = list(telem.groupby('deployment'))
    for grp in grps:
        print('# -- Processing telemetered deployment %s' % grp[0])
        deployments.append(suna_datalogger(grp[1], burst=True))
    deployments = [i for i in deployments if i]
    telem = xr.concat(deployments, 'time')

    print('##### Downloading the recovered_host NUTNR data for %s #####' % site)
    rhost = load_gc_thredds(site, node, sensor, 'recovered_host', stream[1], tag)
    deployments = []
    print('# -- Group the data by deployment and process the data')
    grps = list(rhost.groupby('deployment'))
    for grp in grps:
        print('# -- Processing recovered_host deployment %s' % grp[0])
        deployments.append(suna_datalogger(grp[1], burst=True))
    deployments = [i for i in deployments if i]
    rhost = xr.concat(deployments, 'time')

    print('##### Downloading the recovered_inst NUTNR data for %s #####' % site)
    rinst = load_gc_thredds(site, node, sensor, 'recovered_inst', stream[2], tag)
    deployments = []
    print('# -- Group the data by deployment and process the data')
    grps = list(rinst.groupby('deployment'))
    for grp in grps:
        print('# -- Processing recovered_inst deployment %s' % grp[0])
        deployments.append(suna_instrument(grp[1], burst=True))
    deployments = [i for i in deployments if i]
    rinst = xr.concat(deployments, 'time')

    if site in ['CE01ISSM', 'CE06ISSM']:
        # merge and resample to a two-hour data record
        merged = combine_datasets(telem, rhost, rinst, 120)
    else:
        # merge and resample to an hourly data record
        merged = combine_datasets(telem, rhost, rinst, 60)

    return merged


def generate_qartod(site, node, sensor, cut_off):
    """
    Load all NUTNR data for a defined reference designator (using the site,
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
        QARTOD climatology range tables.
    """
    # load the combined data for the different sources of NUTNR data
    data = combine_delivery_methods(site, node, sensor)

    # create boolean arrays of the data marked as "fail" by the quality checks and generate initial
    # HITL annotations that can be combined with system annotations to create a cleaned up data set
    # prior to calculating the QARTOD test values
    fail = data.nitrate_sensor_quality_flag.where(data.nitrate_sensor_quality_flag > 3).notnull()
    blocks = identify_blocks(fail, [24, 48])
    hitl = create_annotations(site, node, sensor, blocks)

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
    data['corrected_nitrate_concentration'][fail] = np.nan
    if 'rollup_annotations_qc_results' in data.variables:
        fail = data.rollup_annotations_qc_results.where(data.rollup_annotations_qc_results > 3).notnull()
        data['corrected_nitrate_concentration'][fail] = np.nan

    if 'corrected_nitrate_concentration_annotations_qc_results' in data.variables:
        fail = data.corrected_nitrate_concentration_annotations_qc_results.where(
            data.corrected_nitrate_concentration_annotations_qc_results > 3).notnull()
        data['corrected_nitrate_concentration'][fail] = np.nan

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
    parameters = ['corrected_nitrate_concentration']
    limits = [-2.0, 3000]

    # create the initial gross range entry
    gr_lookup = process_gross_range(data, parameters, limits, site=site,
                                    node=node, sensor=sensor, stream='suna')

    if node in ['RID16', 'RID26']:
        # replicate it twice for the different streams
        gr_lookup = pd.concat([gr_lookup] * 2, ignore_index=True)

        # re-work the gross range entries for the different streams
        gr_lookup['stream'][0] = 'suna_dcl_recovered'
        gr_lookup['stream'][1] = 'suna_instrument_recovered'
    else:
        gr_lookup['stream'] = 'nutnr_j_cspp_instrument_recovered'

    # add the stream name and the source comment
    gr_lookup['notes'] = ('User range based on data collected through {}.'.format(src_date))

    # create and format the climatology lookups and tables for the data
    clm_lookup, clm_table = process_climatology(data, parameters, limits, site=site, node=node,
                                                sensor=sensor, stream='suna')

    if node in ['RID16', 'RID26']:
        # replicate it twice for the different streams
        clm_lookup = pd.concat([clm_lookup] * 2, ignore_index=True)

        # re-work the climatology entries for the different streams
        clm_lookup['stream'][0] = 'suna_dcl_recovered'
        clm_lookup['stream'][1] = 'suna_instrument_recovered'
    else:
        clm_lookup['stream'] = 'nutnr_j_cspp_instrument_recovered'


    return annotations, gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the NUTNR data from the Gold Copy THREDDS server and create the
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
    parameters = ['bback', 'estimated_chlorophyll', 'fluorometric_cdom']
    for i in range(len(parameters)):
        tbl = '-'.join([site, node, sensor, parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
