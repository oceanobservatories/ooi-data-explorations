#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the PCO2W data from the uncabled surface moorings, re-processing
    the data using quality flags to identify and remove bad data
"""
import numpy as np
import os
import pandas as pd
import xarray as xr

from cgsn_processing.process.finding_calibrations import find_calibration
from cgsn_processing.process.proc_pco2w import Calibrations

from ooi_data_explorations.common import ENCODINGS, get_annotations, get_vocabulary, get_deployment_dates, \
    load_gc_thredds, add_annotation_qc_flags, update_dataset
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_pco2w import pco2w_datalogger, pco2w_instrument
from ooi_data_explorations.qartod.qc_processing import identify_blocks, create_annotations, inputs


def combine_delivery_methods(site, node, sensor):
    """
    Takes the downloaded data from each of the three data delivery methods and
    combines them into a single, merged xarray data set.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return merged:
    """
    # download the telemetered data and re-process it to create a more useful and coherent data set
    tag = '^(?!.*blank).*PCO2W.*nc$'
    telem = load_gc_thredds(site, node, sensor, 'telemetered', 'pco2w_abc_dcl_instrument', tag)
    print('# Group the telemetered data by deployment and process the data')
    deployments = []
    grps = list(telem.groupby('deployment'))
    for grp in grps:
        print('# -- Group and process telemetered deployment %d' % grp[0])
        data = pco2w_datalogger(grp[1])
        start, end = get_deployment_dates(site, node, sensor, grp[0])
        data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
        deployments.append(data)

    # create the final telemetered data set
    telem = xr.concat(deployments, 'time')
    telem = apply_quality_flags(telem, site, node, sensor)

    # download the recovered host data and re-process it to create a more useful and coherent data set
    rhost = load_gc_thredds(site, node, sensor, 'recovered_host', 'pco2w_abc_dcl_instrument_recovered', tag)
    print('# Group the recovered host data by deployment and process the data')
    deployments = []
    grps = list(rhost.groupby('deployment'))
    for grp in grps:
        print('# -- Group and process recovered_host deployment %d' % grp[0])
        data = pco2w_datalogger(grp[1])
        start, end = get_deployment_dates(site, node, sensor, grp[0])
        data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
        deployments.append(data)

    # create the final recovered_host data set
    rhost = xr.concat(deployments, 'time')
    rhost = apply_quality_flags(rhost, site, node, sensor)

    # download the recovered instrument data and re-process it to create a more useful and coherent data set
    rinst = load_gc_thredds(site, node, sensor, 'recovered_inst', 'pco2w_abc_instrument', tag)
    print('# Group the recovered instrument data by deployment and process the data')
    deployments = []
    grps = list(rinst.groupby('deployment'))
    for grp in grps:
        print('# -- Group and process recovered_inst deployment %d' % grp[0])
        data = pco2w_instrument(grp[1])
        start, end = get_deployment_dates(site, node, sensor, grp[0])
        data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
        deployments.append(data)

    # create the final recovered_inst data set
    rinst = xr.concat(deployments, 'time')
    rinst = apply_quality_flags(rinst, site, node, sensor)

    # combine the three datasets into a single, merged time series resampled to a 3-hour interval time series
    merged = combine_datasets(telem, rhost, rinst, 180)
    return merged


def apply_quality_flags(data, site, node, sensor):
    """
    Use the quality flags, based on criteria set by the vendor for the raw
    data, to identify and remove bad pCO2 data (by replacing with a NaN) prior
    to any subsequent calculations. Additionally, loads any system annotations
    for the instrument and uses time ranges marked as fail by the operators
    to exclude bad data.

    :param data: xarray dataset with vendor-based quality flags added
    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return: xarray dataset with the quality flags used to NaN the seawater pCO2
        values that failed the quality tests.
    """
    # pull in the calibration coefficients from the GitHub Asset Management repo
    cal = Calibrations(None)
    csv_url = find_calibration('PCO2W', data.SerialNumber, (data.time.values.astype('int64') * 10 ** -9)[0])
    cal.read_csv(csv_url)

    # NaN calculated pCO2 values that exceed the cal range by more than 25%
    lower = cal.coeffs['cal_range'][0] - (cal.coeffs['cal_range'][0] * 0.25)
    upper = cal.coeffs['cal_range'][1] + (cal.coeffs['cal_range'][1] * 0.25)
    m = (data.pco2_seawater < lower) | (data.pco2_seawater > upper)
    data.pco2_seawater[m] = np.NaN
    data.pco2_seawater_quality_flag[m] = 4

    # create a boolean array of the data marked as "fail" by the pCO2 quality checks and generate initial
    # HITL annotations that can be combined with system annotations and pCO2 quality checks to create
    # a cleaned up data set prior to calculating the QARTOD test values
    fail = data.pco2_seawater_quality_flag.where(data.pco2_seawater_quality_flag == 4).notnull()
    blocks = identify_blocks(fail, [18, 36])
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

    # create a roll-up annotation flag
    data = add_annotation_qc_flags(data, annotations)

    # clean-up the data, NaN-ing values that fail the pCO2 quality checks or were marked as fail in the annotations
    m = (data.pco2_seawater_quality_flag >= 3) | (data.rollup_annotations_qc_results == 4)
    data.pco2_seawater[m] = np.nan
    return data


def combine_mooring_pco2w(site, node, sensor):
    """
    Download all the pCO2 data for a site and combine it into a 3-hour, median
    averaged dataset with the quality flags used to remove bad data points
    prior to the averaging

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return: xarray dataset with all the pCO2 data for a particular mooring
        combined into a single 3-hour, median averaged data set with the bad
        seawater pCO2 values NaN'ed out.
    """
    # gather and combine all the pCO2 data for a particular site, node sensor combination
    data = combine_delivery_methods(site, node, sensor)

    # update the dataset prior to exporting as a NetCDF file
    vocab = get_vocabulary(site, node, sensor)[0]
    data = update_dataset(data, vocab['maxdepth'])
    return data


def main(argv=None):
    # set up the input arguments
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor

    # create the output directory, if needed
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/adhoc/pco2w')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # gather the data ...
    data = combine_mooring_pco2w(site, node, sensor)

    # ... and save the results to disk
    out_file = ('%s-%s-%s.combined_seawater_pco2.nc' % (site, node, sensor))
    nc_out = os.path.join(out_path, out_file)
    data.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
