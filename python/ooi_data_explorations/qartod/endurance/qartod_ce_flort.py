#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the FLORT data from the uncabled, Coastal Endurance Surface
    Moorings and processes the data to generate QARTOD Gross Range and
    Climatology test limits
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz

from ooi_data_explorations.common import get_annotations, load_gc_thredds, add_annotation_qc_flags
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_flort import flort_datalogger, flort_instrument, flort_cspp, flort_wfp
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology, inputs


def clean_overlaps(ds):
    """
    Clean up the overlapping datasets to make sure we monotonic time records

    :param ds: dataset to sort and reindex
    :return ds: sorted dataset
    """
    ds = ds.sortby('time')
    _, index = np.unique(ds['time'], return_index=True)
    ds = ds.isel(time=index)
    return ds


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

    # download the telemetered data, if available
    telem = load_gc_thredds(site, node, sensor, 'telemetered', stream, tag)

    if node in ['SP001', 'WFP01']:
        # this FLORT is part of a CSPP or WFP and includes telemetered and recovered host data
        if node == 'SP001':
            telem = None  # don't use the telemetered CSPP data
            rhost = load_gc_thredds(site, node, sensor, 'recovered_cspp', stream, tag)
            rhost = flort_cspp(rhost)
        else:
            telem = flort_wfp(telem)
            rhost = load_gc_thredds(site, node, sensor, 'telemetered_wfp', stream, tag)
            rhost = flort_wfp(rhost)

        # merge, but do not resample the time records.
        merged = combine_datasets(telem, rhost, None, None)
    elif node == 'SBD17':
        # this FLORT is mounted on the buoy of the Inshore moorings and includes all three types of data
        telem = flort_datalogger(telem, False)
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host', stream, tag)
        rhost = flort_datalogger(rhost, False)
        rinst = load_gc_thredds(site, node, sensor, 'recovered_inst', stream, tag)
        rinst = flort_instrument(rinst)

        # merge and resample to a 3 hour data record
        merged = combine_datasets(telem, rhost, rinst, 180)
    else:
        # this FLORT is on one of the NSIFs and includes the telemetered and recovered_host data
        telem = clean_overlaps(telem)
        telem = flort_datalogger(telem, True)
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host', stream, tag)
        rhost = clean_overlaps(rhost)
        rhost = flort_datalogger(rhost, True)

        if node == 'RID16':
            # this is an Inshore NSIF, merge and resample to a 3 hour record
            merged = combine_datasets(telem, rhost, None, 180)
        else:
            # this is an coastal NSIF, merge and resample to an hourly record
            merged = combine_datasets(telem, rhost, None, 60)

    return merged


def generate_qartod(site, node, sensor, cut_off):
    """
    Load all of the FLORT data for a defined reference designator (using the
    site, node and sensor names to construct the reference designator) and
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

    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

    # create an annotation-based quality flag
    data = add_annotation_qc_flags(data, annotations)

    # use the existing qc results to NaN values that have failed the different tests
    'bback_qc_results',
    'beta_700_qc_results',
    'estimated_chlorophyll_qc_results',

    # clean-up the data, removing values that were marked as fail in the annotations

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

    # set the parameters and the pressure limits
    parameters = ['seawater_temperature', 'abs_seafloor_pressure',
                  'seawater_temperature', 'abs_seafloor_pressure',
                  'presf_tide_temperature', 'presf_tide_pressure']
    if site in ['CE01ISSM', 'CE06ISSM']:
        plimit = [0, 70]    # 100 psia pressure sensor
    elif site == 'CE07SHSM':
        plimit = [0, 207]   # 300 psia pressure sensor
    else:
        plimit = [0, 689]   # 1000 psia pressure sensor

    limits = [[-5, 35], plimit, [-5, 35], plimit, [-5, 35], plimit]

    # create the initial gross range entry
    gr_lookup = process_gross_range(data, parameters, limits, site=site, node=node, sensor=sensor)

    # re-work gross entry for the different streams
    gr_lookup['stream'][0] = 'presf_abc_dcl_tide_measurement'
    gr_lookup['source'] = ('Sensor min/max based on the vendor sensor specifications. '
                           'The user min/max is the historical mean of all data collected '
                           'up to {} +/- 3 standard deviations.'.format(src_date))

    # create and format the climatology lookups and tables for the data
    clm_lookup, clm_table = process_climatology(data, parameters, limits, site=site, node=node, sensor=sensor)

    # re-work climatology entries for the different streams
    clm_lookup['stream'][0] = 'presf_abc_dcl_tide_measurement'
    clm_lookup['stream'][1] = 'presf_abc_dcl_tide_measurement'
    clm_lookup['stream'][2] = 'presf_abc_dcl_tide_measurement_recovered'
    clm_lookup['stream'][3] = 'presf_abc_dcl_tide_measurement_recovered'
    clm_lookup['stream'][4] = 'presf_abc_tide_measurement_recovered'
    clm_lookup['stream'][5] = 'presf_abc_tide_measurement_recovered'

    return annotations, gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the PRESF data from the Gold Copy THREDDS server and create the
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
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/presf')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # save the annotations to a csv file for further processing
    csv_names = ['id', 'subsite', 'node', 'sensor', 'method', 'stream', 'parameters',
                 'beginDate', 'endDate', 'exclusionFlag', 'qcFlag', 'source', 'annotation']
    anno_csv = '-'.join([site, node, sensor]) + '.quality_annotations.csv'
    annotations.to_csv(os.path.join(out_path, anno_csv), index=False, columns=csv_names)

    # save the gross range values to a csv for further processing
    csv_names = ['subsite', 'node', 'sensor', 'stream', 'parameter', 'qcConfig', 'source']
    gr_csv = '-'.join([site, node, sensor]) + '.gross_range.csv'
    gr_lookup.to_csv(os.path.join(out_path, gr_csv), index=False, columns=csv_names)

    # save the climatology values and table to a csv for further processing
    csv_names = ['subsite', 'node', 'sensor', 'stream', 'parameters', 'climatologyTable', 'source']
    clm_csv = '-'.join([site, node, sensor]) + '.climatology.csv'
    clm_lookup.to_csv(os.path.join(out_path, clm_csv), index=False, columns=csv_names)
    parameters = ['seawater_temperature', 'abs_seafloor_pressure',
                  'seawater_temperature', 'abs_seafloor_pressure',
                  'presf_tide_temperature', 'presf_tide_pressure']
    for i in range(len(parameters)):
        tbl = '-'.join([site, node, sensor, parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
