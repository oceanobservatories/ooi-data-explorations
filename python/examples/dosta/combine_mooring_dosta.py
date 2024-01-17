#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the DOSTA data from the uncabled surface moorings, re-processing
    the data using quality flags to identify and remove bad data
"""
import numpy as np
import os
import pandas as pd
import xarray as xr

from ooi_data_explorations.common import ENCODINGS, get_annotations, get_vocabulary, get_deployment_dates, \
    load_gc_thredds, add_annotation_qc_flags, update_dataset
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_dosta import dosta_datalogger, dosta_ctdbp_datalogger, \
    dosta_ctdbp_instrument
from ooi_data_explorations.qartod.qc_processing import inputs


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
    tag = '.*DOSTA.*nc$'
    if node == 'RID27':
        telem = load_gc_thredds(site, node, sensor, 'telemetered', 'dosta_abcdjm_dcl_instrument', tag)
        print('# Group the telemetered data by deployment and process the data')
        deployments = []
        grps = list(telem.groupby('deployment'))
        for grp in grps:
            print('# -- Group and process telemetered deployment %d' % grp[0])
            data = dosta_datalogger(grp[1], 'True')
            start, end = get_deployment_dates(site, node, sensor, grp[0])
            data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
            deployments.append(data)

        # create the final telemetered data set
        telem = xr.concat(deployments, 'time')
        telem = apply_quality_flags(telem, site, node, sensor)

        # download the recovered host data and re-process it to create a more useful and coherent data set
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host', 'dosta_abcdjm_dcl_instrument_recovered', tag)
        print('# Group the recovered host data by deployment and process the data')
        deployments = []
        grps = list(rhost.groupby('deployment'))
        for grp in grps:
            print('# -- Group and process recovered_host deployment %d' % grp[0])
            data = dosta_datalogger(grp[1], True)
            start, end = get_deployment_dates(site, node, sensor, grp[0])
            data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
            deployments.append(data)

        # create the final recovered_host data set
        rhost = xr.concat(deployments, 'time')
        rhost = apply_quality_flags(rhost, site, node, sensor)

        # no recovered instrument data for this sensor
        rinst = None
    else:
        telem = load_gc_thredds(site, node, sensor, 'telemetered', 'dosta_abcdjm_ctdbp_dcl_instrument', tag)
        print('# Group the telemetered data by deployment and process the data')
        deployments = []
        grps = list(telem.groupby('deployment'))
        for grp in grps:
            print('# -- Group and process telemetered deployment %d' % grp[0])
            data = dosta_ctdbp_datalogger(grp[1])
            start, end = get_deployment_dates(site, node, sensor, grp[0])
            data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
            deployments.append(data)

        # create the final telemetered data set
        telem = xr.concat(deployments, 'time')
        telem = apply_quality_flags(telem, site, node, sensor)

        # download the recovered host data and re-process it to create a more useful and coherent data set
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host', 'dosta_abcdjm_ctdbp_dcl_instrument_recovered', tag)
        print('# Group the recovered host data by deployment and process the data')
        deployments = []
        grps = list(rhost.groupby('deployment'))
        for grp in grps:
            print('# -- Group and process recovered_host deployment %d' % grp[0])
            data = dosta_ctdbp_datalogger(grp[1])
            start, end = get_deployment_dates(site, node, sensor, grp[0])
            data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
            deployments.append(data)

        # create the final recovered_host data set
        rhost = xr.concat(deployments, 'time')
        rhost = apply_quality_flags(rhost, site, node, sensor)

        # download the recovered instrument data and re-process it to create a more useful and coherent data set
        rinst = load_gc_thredds(site, node, sensor, 'recovered_inst', 'dosta_abcdjm_ctdbp_instrument_recovered', tag)
        print('# Group the recovered instrument data by deployment and process the data')
        deployments = []
        grps = list(rinst.groupby('deployment'))
        for grp in grps:
            print('# -- Group and process recovered_inst deployment %d' % grp[0])
            data = dosta_ctdbp_instrument(grp[1])
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
    Load and system annotations for the instrument and use time ranges marked
    as fail by the operators to exclude bad data.

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
    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

    # create a roll-up annotation flag
    data = add_annotation_qc_flags(data, annotations)

    # clean-up the data, NaN-ing values that were marked as fail in the annotations
    m = (data.rollup_annotations_qc_results == 4)
    data.oxygen_concentration[m] = np.nan
    data.oxygen_concentration_corrected[m] = np.nan

    # catch any system fill values (-999999999) and NaN them out
    m = (data.oxygen_concentration_corrected < -10)
    data.oxygen_concentration[m] = np.nan
    data.oxygen_concentration_corrected[m] = np.nan

    return data


def combine_mooring_dosta(site, node, sensor):
    """
    Download all the dissolved oxygen data for a site and combine it into a
    3-hour, median averaged dataset with the quality flags used to remove bad
    data points prior to the averaging

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return: xarray dataset with all the DO data for a particular mooring
        combined into a single 3-hour, median averaged data set with the bad
        values NaN'ed out.
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
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/adhoc/dosta')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # gather the data ...
    data = combine_mooring_dosta(site, node, sensor)

    # ... and save the results to disk
    out_file = ('%s-%s-%s.combined_dissolved_oxygen.nc' % (site, node, sensor))
    nc_out = os.path.join(out_path, out_file)
    data.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
