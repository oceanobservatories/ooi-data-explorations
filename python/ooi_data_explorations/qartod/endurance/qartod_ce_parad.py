#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the PARAD data from the uncabled, Coastal Endurance Profiler
    Mooring and Coastal Surface-Piercing Profilers (CSPP) and process the data
    to generate QARTOD Gross Range and Climatology test limits
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz
import xarray as xr

from ooi_data_explorations.common import get_annotations, get_vocabulary, load_gc_thredds, add_annotation_qc_flags
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_parad import parad_cspp, parad_wfp
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology, woa_standard_bins, \
    inputs, ANNO_HEADER, CLM_HEADER, GR_HEADER
from ooi_data_explorations.qartod.endurance.qartod_ce_gliders import collect_glider

def combine_delivery_methods(site, node, sensor):
    """
    Takes the downloaded data from the different data delivery methods for the
    profiling Photosynthetically Active Radiation (PARAD) sensors, and combines
    them into a single, merged xarray data set.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return merged: the merged PARAD dataset
    """
    # set the tag constant
    tag = '.*PARAD.*\\.nc$'

    if node == 'SP001':
        if site in ['CE01ISSP', 'CE06ISSP']:
            # this PARAD is part of a CSPP and includes recovered data only
            print('##### Downloading the recovered_cspp PARAD data for %s #####' % site)
            rhost = load_gc_thredds(site, node, sensor, 'recovered_cspp', 'parad_j_cspp_instrument_recovered', tag)
            deployments = []
            print('# -- Group the data by deployment and process the data')
            grps = list(rhost.groupby('deployment'))
            for grp in grps:
                print('# -- Processing recovered_host deployment %s' % grp[0])
                deployments.append(parad_cspp(grp[1]))
            deployments = [i for i in deployments if i]
            rhost = xr.concat(deployments, 'time')

            # merge, but do not resample the time records.
            merged = combine_datasets(None, rhost, None, None)
        else:
            # use the glider data records to develop the QARTOD test limits
            if site == 'CE02SHSP':
                merged = collect_glider(44.639, -124.304)
            else:
                merged = collect_glider(46.986, -124.566)

    else:
        # this PARAD is part of a WFP and includes telemetered and recovered data
        print('##### Downloading the telemetered PARAD data for %s #####' % site)
        telem = load_gc_thredds(site, node, sensor, 'telemetered', 'parad_k__stc_imodem_instrument', tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(telem.groupby('deployment'))
        for grp in grps:
            print('# -- Processing telemetered deployment %s' % grp[0])
            deployments.append(parad_wfp(grp[1]))
        deployments = [i for i in deployments if i]
        telem = xr.concat(deployments, 'time')

        print('##### Downloading the recovered_wfp PARAD data for %s #####' % site)
        rhost = load_gc_thredds(site, node, sensor, 'recovered_wfp', 'parad_k__stc_imodem_instrument_recovered', tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(rhost.groupby('deployment'))
        for grp in grps:
            print('# -- Processing recovered_host deployment %s' % grp[0])
            deployments.append(parad_wfp(grp[1]))
        deployments = [i for i in deployments if i]
        rhost = xr.concat(deployments, 'time')

        # merge, but do not resample the time records.
        merged = combine_datasets(telem, rhost, None, None)

    return merged


def generate_qartod(site, node, sensor, cut_off):
    """
    Load all PARAD data for a defined reference designator (using the site,
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
    # load the combined data for the different sources of PARAD data
    data = combine_delivery_methods(site, node, sensor)

    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty and site not in ['CE02SHSP', 'CE07SHSP']:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

        # create an annotation-based quality flag
        data = add_annotation_qc_flags(data, annotations)

    # clean-up the data, NaN-ing values that were marked as fail in the QC checks and then removing all
    # records where the rollup annotation (every parameter fails) was set to fail.
    if 'par_qc_summary_flag' in data.variables:
        m = data.par_qc_summary_flag == 4
        data['par'][m] = np.nan

    if 'rollup_annotations_qc_results' in data.variables:
        data = data.where(data.rollup_annotations_qc_results != 4, drop=True)

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

    _, index = np.unique(data['time'], return_index=True)
    data = data.isel(time=index)
    data = data.sel(time=slice('2014-01-01T00:00:00', end_date))

    # set the parameters and the gross range limits
    parameters = ['par']
    limits = [0, 5000]

    # create the initial gross range entry
    quantile = data['depth'].quantile([0.01, 0.99], dim='time', interpolation='linear')
    sub = data.where(data.depth <= quantile[0], drop=True)  # limit gross range estimation to near-surface values
    gr_lookup = process_gross_range(sub, parameters, limits, site=site,
                                    node=node, sensor=sensor, stream='parad_replace_me')

    # add the stream name and the source comment
    gr_lookup['notes'] = ('User range based on data collected through {}.'.format(src_date))

    # based on node, determine if we need a depth based climatology
    depth_bins = np.array([])
    if node == 'WFP01':
        vocab = get_vocabulary(site, node, sensor)[0]
        max_depth = vocab['maxdepth']
        depth_bins = woa_standard_bins()
        m = depth_bins[:, 1] <= max_depth
        depth_bins = depth_bins[m, :]
        depth_bins = depth_bins[6:, :]  # profiler doesn't go above 30 m

    if node == 'SP001' and site in ['CE02SHSP', 'CE07SHSP']:
        vocab = get_vocabulary(site, node, sensor)[0]
        max_depth = vocab['maxdepth']
        depth_bins = woa_standard_bins()
        m = depth_bins[:, 1] <= max_depth
        depth_bins = depth_bins[m, :]

    # create and format the climatology lookups and tables for the data
    clm_lookup, clm_table = process_climatology(data, parameters, limits, depth_bins=depth_bins,
                                                site=site, node=node, sensor=sensor,
                                                stream='parad_replace_me')

    return annotations, gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the PARAD data from the Gold Copy THREDDS server and create the
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
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/parad')
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
    parameters = ['par']
    for i in range(len(parameters)):
        tbl = '-'.join([site, node, sensor, parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
