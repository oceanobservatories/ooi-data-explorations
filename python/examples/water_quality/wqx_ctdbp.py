#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the CTDBP data from the uncabled, Coastal Endurance Surface
    Moorings and processes the data to generate a data file for the Oregon DEQ
    water quality project (https://www.oregon.gov/deq/wq/Pages/index.aspx)
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd

from ooi_data_explorations.common import get_annotations, load_gc_thredds, add_annotation_qc_flags
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.qartod.qc_processing import ANNO_HEADER, inputs
from ooi_data_explorations.uncabled.process_ctdbp import ctdbp_datalogger, ctdbp_instrument


def combine_delivery_methods(site, node, sensor, resample=180):
    """
    Takes the downloaded data from each of the three data delivery methods for
    the uncabled CTD (CTDBP), and combines each of them into a single, merged
    xarray data set.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :param resample: integer value representing the time interval (in minutes)
        to resample the data to (default is 180 minutes)
    :return merged: the merged CTDBP data stream resampled to a user-defined
        time interval
    """
    # set the parameter names of interest
    parameters = ['sea_water_electrical_conductivity', 'sea_water_temperature', 'sea_water_pressure',
                  'sea_water_practical_salinity']

    # download the telemetered data and re-process it to create a more useful and coherent data set
    tag = '.*CTDBP.*\\.nc$'
    telem = load_gc_thredds(site, node, sensor, 'telemetered', 'ctdbp_cdef_dcl_instrument', tag)
    telem = ctdbp_datalogger(telem)

    # download the recovered host data and re-process it to create a more useful and coherent data set
    rhost = load_gc_thredds(site, node, sensor, 'recovered_host', 'ctdbp_cdef_dcl_instrument_recovered', tag)
    rhost = ctdbp_datalogger(rhost)

    # download the recovered instrument data and re-process it to create a more useful and coherent data set
    rinst = load_gc_thredds(site, node, sensor, 'recovered_inst', 'ctdbp_cdef_instrument_recovered', tag)
    rinst = ctdbp_instrument(rinst)

    # use the QARTOD results variables to remove data that failed the QARTOD tests
    for param in parameters:
        m = telem[param + '_qartod_results'] >= 4
        telem[param][m] = np.nan
        m = rhost[param + '_qartod_results'] >= 4
        rhost[param][m] = np.nan
        m = rinst[param + '_qartod_results'] >= 4
        rinst[param][m] = np.nan

    # combine the three datasets into a single, merged time series resampled to a user defined time interval (default
    # is 3 hours)
    merged = combine_datasets(telem, rhost, rinst, resample_time=resample)
    return merged


def generate_deq(site, node, sensor, resample, cut_off):
    """
    Load the CTD data for a defined reference designator (using the site, node
    and sensor names to construct the reference designator) collected via the
    telemetered, recovered host and instrument methods and combine them into a
    single data set that can be used to create a data file for the DEQ water
    quality project.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :param resample: integer value representing the time interval (in
        minutes) to resample the data to. Default is 180 minutes (3 hours).
    :param cut_off: string formatted dates to use as cut-offs for the data
        (e.g. '2018-01-01')
    :return data: xarray data set containing the combined data from the three
        data delivery methods for the CTD instrument
    """
    # load the combined telemetered, recovered_host and recovered_inst data
    data = combine_delivery_methods(site, node, sensor, resample)

    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

    fail = annotations[annotations.qcFlag == 'fail']
    if not fail.empty:
        data = add_annotation_qc_flags(data, fail.reset_index())
        data = data.where(data.rollup_annotations_qc_results < 3, drop=True)

    # slice the data to the user defined cut-off dates
    strt = parser.parse(cut_off[0])
    strt_date = strt.strftime('%Y-%m-%dT%H:%M:%S')
    end = parser.parse(cut_off[1])
    end_date = end.strftime('%Y-%m-%dT%H:%M:%S')
    data = data.sel(time=slice(strt_date, end_date))

    return annotations, data


def main(argv=None):
    """
    Download the CTDBP data from the Gold Copy THREDDS server and create the
    QARTOD gross range and climatology test lookup tables.
    """
    # set up the input arguments
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    cut_off = args.cut_off

    # collect the data, resampling to 60 minutes and remove data outside the user defined cut-off dates
    annotations, data = generate_deq(site, node, sensor, 60, cut_off)

    # save the downloaded annotations
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/wqx/ctdbp')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # save the annotations to a csv file for later review
    anno_csv = '-'.join([site, node, sensor]) + '.quality_annotations.csv'
    annotations.to_csv(os.path.join(out_path, anno_csv), index=False, columns=ANNO_HEADER)

    # save the data to a netCDF file for subsequent processing
    ENCODING = {'time': {'units': 'seconds since 1900-01-01T00:00:00.000Z', 'calendar': 'gregorian'}}
    data.to_netcdf(os.path.join(out_path, '-'.join([site, node, sensor]) + '.nc'), encoding=ENCODING)


if __name__ == '__main__':
    main()
