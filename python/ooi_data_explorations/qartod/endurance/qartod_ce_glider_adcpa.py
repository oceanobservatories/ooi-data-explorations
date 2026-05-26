#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the ADCPA data from the OOI Gold Copy THREDDS server for the
    OOI Coastal Endurance gliders and process the data to generate QARTOD
    Gross Range test limits.

Note: climatology tests are deferred -- the combination of glider platform
motion and depth-variable profiling makes the bin-based depth climatology
approach used for moored ADCPs unsuitable without further methodology work.
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz

from ooi_data_explorations.common import get_annotations, load_gc_thredds, add_annotation_qc_flags
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.qartod.qc_processing import process_gross_range, inputs, ANNO_HEADER, GR_HEADER


def combine_delivery_methods(site: str, node: str, sensor: str):
    """
    Download the ADCPA velocity data from the Gold Copy THREDDS server.

    Glider ADCPs have only recovered_host data; no telemetered or
    recovered_inst streams exist for this instrument class.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and
        fourth part of the reference designator
    :return: merged ADCPA dataset
    """
    tag = '.*ADCP.*\\.nc$'
    print('##### Downloading the recovered_host ADCPA data for %s #####' % site)
    rhost = load_gc_thredds(site, node, sensor, 'recovered_host', 'adcp_velocity_earth', tag)
    merged = combine_datasets(None, rhost, None, None)
    return merged


def generate_qartod(site: str, node: str, sensor: str, cut_off: str) -> tuple:
    """
    Load all ADCPA data for a defined reference designator and combine
    into a single data set from which QARTOD gross range test limits can
    be calculated.

    Climatology tests are not generated for glider ADCPA data; see the
    module docstring for rationale.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and
        fourth part of the reference designator
    :param cut_off: string formatted date to use as cut-off for data to
        add to QARTOD test sets
    :return annotations: DataFrame of system annotations for the sensor
    :return gr_lookup: CSV formatted strings for the QARTOD gross range
        lookup tables
    """
    data = combine_delivery_methods(site, node, sensor)

    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(
            annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(
            annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

        data = add_annotation_qc_flags(data, annotations)

    if 'rollup_annotations_qc_results' in data.variables:
        data = data.where(data.rollup_annotations_qc_results != 4, drop=True)

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

    # set the parameters and sensor range limits
    parameters = ['eastward_seawater_velocity', 'northward_seawater_velocity']
    limits = [[-5, 5], [-5, 5]]

    gr_lookup = process_gross_range(data, parameters, limits, site=site, node=node,
                                    sensor=sensor, stream='adcp_velocity_earth', extended=True)
    gr_lookup['notes'] = ('User range based on data collected through {}.'.format(src_date))

    return annotations, gr_lookup


def main(argv=None):
    """
    Download the ADCPA data from the Gold Copy THREDDS server and create
    the QARTOD gross range test lookup tables.
    """
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    cut_off = args.cut_off

    annotations, gr_lookup = generate_qartod(site, node, sensor, cut_off)

    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/adcpa')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    anno_csv = '-'.join([site, node, sensor]) + '.quality_annotations.csv'
    annotations.to_csv(os.path.join(out_path, anno_csv), index=False, columns=ANNO_HEADER)

    gr_csv = '-'.join([site, node, sensor]) + '.gross_range.csv'
    gr_lookup.to_csv(os.path.join(out_path, gr_csv), index=False, columns=GR_HEADER)


if __name__ == '__main__':
    main()
