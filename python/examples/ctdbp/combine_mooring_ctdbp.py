#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the CTDBP data from the uncabled surface moorings, re-processing
    the data using quality flags to identify and remove bad data
"""
import numpy as np
import os
import pandas as pd
import xarray as xr

from ooi_data_explorations.common import ENCODINGS, get_annotations, get_vocabulary, get_deployment_dates, \
    load_gc_thredds, add_annotation_qc_flags, update_dataset
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_ctdbp import ctdbp_datalogger, ctdbp_instrument
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
    # determine the minimum depth for the sensor and use it to remove data shallower than that depth, helps to remove
    # data collected during deployment and recovery when the instrument is either out of the water or still descending
    # or ascending through the water column
    if node == 'SBD17':
        depth_cutoff = 0.0  # surface buoy nominal depth is 1.5 m
    elif node in ['RID16', 'RID27']:
        depth_cutoff = 3.5  # NSIF nominal depth is 7 m
    else:
        vocab = get_vocabulary(site, node, sensor)
        depth_cutoff = vocab[0]['mindepth'] / 1.33

    # download the telemetered data and re-process it to create a more useful and coherent data set
    tag = '.*CTDBP.*nc$'
    telem = load_gc_thredds(site, node, sensor, 'telemetered', 'ctdbp_cdef_dcl_instrument', tag)
    print('# Group the telemetered data by deployment and process the data')
    deployments = []
    grps = list(telem.groupby('deployment'))
    for grp in grps:
        print('# -- Group and process telemetered deployment %d' % grp[0])
        data = ctdbp_datalogger(grp[1])
        start, end = get_deployment_dates(site, node, sensor, grp[0])
        data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
        data = data.where(data.sea_water_pressure > depth_cutoff, drop=True)  # remove data above the minimum depth
        deployments.append(data)

    # create the final telemetered data set
    telem = xr.concat(deployments, 'time')
    telem = apply_quality_flags(telem, site, node, sensor)

    # download the recovered host data and re-process it to create a more useful and coherent data set
    rhost = load_gc_thredds(site, node, sensor, 'recovered_host', 'ctdbp_cdef_dcl_instrument_recovered', tag)
    print('# Group the recovered host data by deployment and process the data')
    deployments = []
    grps = list(rhost.groupby('deployment'))
    for grp in grps:
        print('# -- Group and process recovered_host deployment %d' % grp[0])
        data = ctdbp_datalogger(grp[1])
        start, end = get_deployment_dates(site, node, sensor, grp[0])
        data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
        data = data.where(data.sea_water_pressure > depth_cutoff, drop=True)  # remove data above the minimum depth
        deployments.append(data)

    # create the final recovered_host data set
    rhost = xr.concat(deployments, 'time')
    rhost = apply_quality_flags(rhost, site, node, sensor)

    # download the recovered instrument data and re-process it to create a more useful and coherent data set
    rinst = load_gc_thredds(site, node, sensor, 'recovered_inst', 'ctdbp_cdef_instrument_recovered', tag)
    print('# Group the recovered instrument data by deployment and process the data')
    deployments = []
    grps = list(rinst.groupby('deployment'))
    for grp in grps:
        print('# -- Group and process recovered_inst deployment %d' % grp[0])
        data = ctdbp_instrument(grp[1])
        start, end = get_deployment_dates(site, node, sensor, grp[0])
        data = data.sel(time=slice(start[:-1], end[:-1]))  # trim the data to the deployment dates
        data = data.where(data.sea_water_pressure > depth_cutoff, drop=True)  # remove data above the minimum depth
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
    data, to identify and remove bad CTD data (by replacing with a NaN) prior
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
    :return: xarray dataset with the quality flags used to NaN the CTD data
        values that failed the quality tests.
    """
    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

    # add roll-up and variable specific annotation flags to the data
    data = add_annotation_qc_flags(data, annotations)

    # clean-up the data, removing data where the roll-up flag is fail, and NaN-ing values that fail the different
    # quality checks (automated and HITL)
    data = data.where(data.rollup_annotations_qc_results != 4)

    qartod = ['sea_water_electrical_conductivity_qartod_results', 'sea_water_temperature_qartod_results',
              'sea_water_pressure_qartod_results', 'sea_water_practical_salinity_qartod_results']
    for q in qartod:
        m = data[q] == 4
        data[q.replace('_qartod_results', '')][m] = np.nan
        # propagate the NaNs to the other variables that are calculated from the variable that failed the QC
        if q in ['sea_water_electrical_conductivity_qartod_results', 'sea_water_temperature_qartod_results',
                 'sea_water_pressure_qartod_results']:
            data['sea_water_practical_salinity'][m] = np.nan
            data['sea_water_density'][m] = np.nan
        else:
            data['sea_water_density'][m] = np.nan

    hitl = ['sea_water_electrical_conductivity_annotations_qc_results', 'sea_water_temperature_annotations_qc_results',
            'sea_water_pressure_annotations_qc_results']
    for h in hitl:
        if h in data.variables:
            m = data[h] == 4
            data[h.replace('_annotations_qc_results', '')][m] = np.nan
            # propagate the NaNs to the other variables that are calculated from the variable that failed the QC
            if h in ['sea_water_electrical_conductivity_annotations_qc_results',
                     'sea_water_temperature_annotations_qc_results', 'sea_water_pressure_annotations_qc_results']:
                data['sea_water_practical_salinity'][m] = np.nan
                data['sea_water_density'][m] = np.nan
            else:
                data['sea_water_density'][m] = np.nan

    return data


def combine_mooring_ctdbp(site, node, sensor):
    """
    Download all the CTD data for a site and combine it into a 3-hour, median
    averaged dataset with the quality flags used to remove bad data points
    prior to the averaging

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return: xarray dataset with all the CTD data for a particular mooring
        combined into a single 3-hour, median averaged data set with the bad
        CTD values NaN'ed out.
    """
    # gather and combine all the CTD data for a particular site, node and sensor combination
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
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/adhoc/ctdbp')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # gather the data ...
    data = combine_mooring_ctdbp(site, node, sensor)

    # ... and save the results to disk
    out_file = ('%s-%s-%s.combined_ctdbp.nc' % (site, node, sensor))
    nc_out = os.path.join(out_path, out_file)
    data.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
