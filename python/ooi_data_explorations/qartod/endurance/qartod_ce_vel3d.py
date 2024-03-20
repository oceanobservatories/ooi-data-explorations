#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the VEL3D data from the uncabled, Coastal Endurance Surface
    Moorings and Profilers and process the data to generate QARTOD Gross Range
    and Climatology test limits
"""
import dateutil.parser as parser
import numpy as np
import os
import pandas as pd
import pytz
import xarray as xr

from ooi_data_explorations.common import get_annotations, get_vocabulary, load_gc_thredds, get_deployment_dates, \
    m2m_request, m2m_collect, add_annotation_qc_flags
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_vel3d import vel3d_datalogger, mmp_aquadopp
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology, \
    woa_standard_bins, inputs, ANNO_HEADER, CLM_HEADER, GR_HEADER


def combine_delivery_methods(site, node, sensor):
    """
    Takes the downloaded data from the different data delivery methods for the
    seven-channel, downwelling spectral irradiance (VEL3D) sensor, and combines
    them into a single, merged xarray data sets.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return merged: the merged and resampled (if appropriate) VEL3D dataset
    """
    # set the stream and tag constants
    if node == 'WFP01':
        tag = '.*VEL3D.*\\.nc$'
        # only download the recovered data for the MMP Aquadopp II, the telemetered data is not that useful
        print('##### Downloading the recovered Aquadopp II data for %s #####' % site)
        rinst = load_gc_thredds(site, node, sensor, 'recovered_wfp', 'vel3d_k_wfp_instrument', tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(rinst.groupby('deployment'))
        for grp in grps:
            print('# -- Processing recovered deployment %s' % grp[0])
            deployments.append(mmp_aquadopp(grp[1], binning=True, bin_size=2.0))

        # concatenate the deployments into a single dataset
        deployments = [i for i in deployments if i]
        merged = xr.concat(deployments, 'time')
        merged = merged.sortby(['deployment', 'profile', 'time'])
        _, index = np.unique(merged['time'], return_index=True)
        merged = merged.isel(time=index)
    else:
        print('##### Downloading the telemetered VEL3D data for %s #####' % site)
        tag = '.*VEL3D.*velocity.*\\.nc$'
        telem = load_gc_thredds(site, node, sensor, 'telemetered', 'vel3d_cd_dcl_velocity_data', tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(telem.groupby('deployment'))
        for grp in grps:
            print('# -- Collecting the telemetered header and system packet data for deployment %s' % grp[0])
            deploy = grp[0]
            tag = '.*deployment%04d.*VEL3D.*system.*\\.nc$' % deploy
            start, stop = get_deployment_dates(site, node, sensor, deploy)
            m = m2m_request(site, node, sensor, 'telemetered', 'vel3d_cd_dcl_system_data', start, stop)
            system = m2m_collect(m, tag)
            tag = '.*deployment%04d.*VEL3D.*header.*\\.nc$' % deploy
            m = m2m_request(site, node, sensor, 'telemetered', 'vel3d_cd_dcl_data_header', start, stop)
            header = m2m_collect(m, tag)
            print('# -- Merging the header, system, and velocity data packets for deployment %s' % deploy)
            deployments.append(vel3d_datalogger(header, system, grp[1], burst=True))

        # concatenate the deployments into a single dataset
        deployments = [i for i in deployments if i]
        telem = xr.concat(deployments, 'time')

        print('##### Downloading the recovered_host VEL3D data for %s #####' % site)
        tag = '.*VEL3D.*velocity.*\\.nc$'
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host', 'vel3d_cd_dcl_velocity_data_recovered', tag)
        deployments = []
        print('# -- Group the data by deployment and process the data')
        grps = list(rhost.groupby('deployment'))
        for grp in grps:
            print('# -- Collecting the recovered_host header and system packet data for deployment %s' % grp[0])
            deploy = grp[0]
            tag = '.*deployment%04d.*VEL3D.*system.*\\.nc$' % deploy
            start, stop = get_deployment_dates(site, node, sensor, deploy)
            m = m2m_request(site, node, sensor, 'recovered_host', 'vel3d_cd_dcl_system_data_recovered', start, stop)
            system = m2m_collect(m, tag)
            tag = '.*deployment%04d.*VEL3D.*header.*\\.nc$' % deploy
            m = m2m_request(site, node, sensor, 'recovered_host', 'vel3d_cd_dcl_data_header_recovered', start, stop)
            header = m2m_collect(m, tag)
            print('# -- Merging the header, system, and velocity data packets for deployment %s' % deploy)
            deployments.append(vel3d_datalogger(header, system, grp[1], burst=True))

        # concatenate the deployments into a single dataset
        deployments = [i for i in deployments if i]
        rhost = xr.concat(deployments, 'time')

        # combine the datasets, leaving them as burst averaged datasets
        merged = combine_datasets(telem, rhost, None, None)

    return merged


def generate_qartod(site, node, sensor, cut_off):
    """
    Load all VEL3D data for a defined reference designator (using the site,
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
    # load the combined data for the different sources of VEL3D data
    data = combine_delivery_methods(site, node, sensor)

    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

    # create an annotation-based quality flag and remove any data that has been flagged as bad
    data = add_annotation_qc_flags(data, annotations)
    if 'rollup_annotations_qc_results' in data.variables:
        data = data.where(data.rollup_annotations_qc_results != 4)

    # if this is an MMP Aquadopp II, rename the QC variable to match the other renamed variables and apply variable
    # specific QC test results to the data, otherwise just apply the variable specific QC test results to the data
    # for the Vector.
    if node == 'WFP01':
        rename = {
            'vel3d_k_upward_velocity_annotations_qc_results': 'relative_velocity_vertical_annotations_qc_results',
            'vel3d_k_eastward_velocity_annotations_qc_results': 'relative_velocity_east_annotations_qc_results',
            'vel3d_k_northward_velocity_annotations_qc_results': 'relative_velocity_north_annotations_qc_results',
            'vel3d_k_data_set_description_annotations_qc_results': 'beam_mapping_annotation_qc_results',
            'vel3d_k_vel0_annotations_qc_results': 'velocity_1_annotations_qc_results',
            'vel3d_k_vel1_annotations_qc_results': 'velocity_2_annotations_qc_results'
        }
        for key in rename.keys():
            if key in data.variables:
                data = data.rename({key: rename.get(key)})

        # apply the variable specific QC test results to the data
        data = data.where(data.aquadopp_sensor_quality_flag != 4)
    else:
        data = data.where(data.vector_sensor_quality_flag != 4)

    # if a cut_off date was used, limit data to all data collected up to the cut_off date.
    # otherwise, set the limit to the range of the downloaded data.
    if cut_off:
        cut = parser.parse(cut_off)
        cut = cut.astimezone(pytz.utc)
        end_date = cut.strftime('%Y-%m-%dT23:59:59.999999')
        src_date = cut.strftime('%Y-%m-%d')
    else:
        cut = parser.parse(data.time_coverage_end)
        cut = cut.astimezone(pytz.utc)
        end_date = cut.strftime('%Y-%m-%dT23:59:59.999999')
        src_date = cut.strftime('%Y-%m-%d')

    data = data.sel(time=slice('2014-01-01T00:00:00', end_date))
    start_date = str(data.time[0].values.min())[:10]

    if node == 'WFP01':
        # we are working with the Aquadopp II on the MMP. Parameters for QARTOD are different for
        # the MMP Aquadopp II than for the Vector.
        parameters = ['sea_water_temperature', 'relative_velocity_east',
                      'relative_velocity_north', 'relative_velocity_vertical']
        limits = [[-400, 4000], [-3.0, 3.0], [-3.0, 3.0], [-1.0, 1.0]]

        # make sure we are using the same scaling as found in the source data
        data['sea_water_pressure'] = data['sea_water_pressure'] / 0.001  # reported in 0.001 dbar
        data['sea_water_temperature'] = data['sea_water_temperature'] / 0.01  # reported in 0.01 degrees Celsius
    else:
        # this is the Vector, set parameters and the gross range limits
        parameters = ['sea_water_pressure', 'sea_water_temperature', 'velocity_east_corrected',
                      'velocity_north_corrected', 'velocity_vertical']
        limits = [[0, 4000000], [-400, 4000], [-2.1, 2.1], [-2.1, 2.1], [-0.6, 0.6]]
        # note, velocity range limits are set based on a software selected nominal range of +/- 1.0 m/s which yields
        # the vertical and horizontal velocity limits of +/- 0.6 m/s and +/- 2.1 m/s shown above.

        # make sure we are using the same scaling as found in the source data
        data['sea_water_pressure'] = data['sea_water_pressure'] / 0.001  # reported in 0.001 dbar
        data['sea_water_temperature'] = data['sea_water_temperature'] / 0.01  # reported in 0.01 degrees Celsius

    # create the initial gross range entry
    gr_lookup = process_gross_range(data, parameters, limits, site=site,
                                    node=node, sensor=sensor, stream='temporary_vel3d_stream_name')

    # add the source comment
    gr_lookup['source'] = ('User Gross Range based on data collected from {} through to {}.'.format(start_date,
                                                                                                    src_date))

    # based on the node, determine if we need a depth based climatology
    depth_bins = np.array([])
    if node == 'WFP01':
        vocab = get_vocabulary(site, node, sensor)[0]
        max_depth = vocab['maxdepth']
        depth_bins = woa_standard_bins()
        m = depth_bins[:, 1] <= max_depth
        depth_bins = depth_bins[m, :]

        # reset parameters, no climatology test for the relative velocity data
        parameters = ['sea_water_temperature']
        limits = [[-400, 4000]]
    else:
        # reset parameters, no climatology test for the pressure data
        parameters = ['sea_water_temperature', 'velocity_east_corrected', 'velocity_north_corrected',
                      'velocity_vertical']
        limits = [[-400, 4000], [-2.1, 2.1], [-2.1, 2.1], [-0.6, 0.6]]

    # create and format the climatology lookups and tables for the data
    clm_lookup, clm_table = process_climatology(data, parameters, limits, depth_bins=depth_bins,
                                                site=site, node=node, sensor=sensor,
                                                stream='temporary_vel3d_stream_name')

    # add the source comment
    clm_lookup['source'] = ('Climatology based on data collected from {} through to {}.'.format(start_date, src_date))

    # return the data, annotations, and the gross range and climatology lookup tables
    return data, annotations, gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the VEL3D data from the Gold Copy THREDDS server and create the
    QARTOD gross range and climatology test lookup tables.
    """
    # set up the input arguments
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    cut_off = args.cut_off

    # create the QARTOD gross range and climatology lookup values and tables
    data, annotations, gr_lookup, clm_lookup, clm_table = generate_qartod(site, node, sensor, cut_off)

    # save the downloaded annotations and qartod lookups and tables
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/vel3d')
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
    if node == 'WFP01':
        parameters = ['sea_water_temperature']
    else:
        # reset parameters, no climatology test for the pressure data
        parameters = ['sea_water_temperature', 'velocity_east_corrected', 'velocity_north_corrected',
                      'velocity_vertical']

    for i in range(len(parameters)):
        tbl = '-'.join([site, node, sensor, parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
