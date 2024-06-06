#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the CTDBP data from the uncabled, Coastal Endurance Surface
    Moorings and processes the data to generate report plots and tables for
    HITL QC review.
"""
import dateutil.parser as parser
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import pandas as pd
import warnings

from matplotlib.backends.backend_pdf import PdfPages
from ooi_data_explorations.common import get_annotations, get_deployment_dates, load_gc_thredds
from ooi_data_explorations.combine_data import combine_datasets
from ooi_data_explorations.uncabled.process_ctdbp import ctdbp_datalogger, ctdbp_instrument
from ooi_data_explorations.qartod.discrete_samples import get_discrete_samples, distance_to_cast
from ooi_data_explorations.qartod.qc_processing import ANNO_HEADER, inputs
from ooi_data_explorations.qartod.reporting import apply_qc_results


def combine_delivery_methods(site, node, sensor, deployment):
    """
    Downloads the  data from each of the three data delivery methods for
    the uncabled CTD (CTDBP) instrument for a specific deployment and
    returns a list object containing the data from each delivery method.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :param deployment: deployment number to download the data for
    :return merged: the three CTDBP data streams as a list object
    """
    # set the regex tag for the data files to download
    tag = 'deployment{:04d}.*CTDBP.*\\.nc$'.format(deployment)
    print('### Downloading data for {}-{}-{}, deployment {:02d} ###'.format(site, node, sensor, deployment))
    warnings.filterwarnings("error")  # set warnings to raise exceptions so we can catch them
    # download the data from the telemetered data stream
    try:
        telem = load_gc_thredds(site, node, sensor, 'telemetered', 'ctdbp_cdef_dcl_instrument', tag)
    except UserWarning as e:
        telem = None
    # download the data from the recovered_host data stream
    try:
        rhost = load_gc_thredds(site, node, sensor, 'recovered_host', 'ctdbp_cdef_dcl_instrument_recovered', tag)
    except UserWarning as e:
        rhost = None
    # download the data from the recovered_inst data stream
    try:
        rinst = load_gc_thredds(site, node, sensor, 'recovered_inst', 'ctdbp_cdef_instrument_recovered', tag)
    except UserWarning as e:
        rinst = None
    # reset the warnings to default behavior
    warnings.resetwarnings()

    # combine the data from the three delivery methods
    return [telem, rhost, rinst]


def plotting(site, node, sensor, deployment, dates, availability, merged, qartod, pre_merged, post_merged, samples):
    """
    Generate plots of the CTDBP data for the telemetered, recovered host and
    recovered instrument data streams, as well as the pre- and post-deployment
    data for the sensor.
    """
    params = {
        'font.size': 10,
        'axes.linewidth': 1.0,
        'axes.titlelocation': 'right',
        'figure.figsize': [17, 11],
        'figure.dpi': 100,
        'xtick.major.size': 4,
        'xtick.major.pad': 4,
        'xtick.major.width': 1.0,
        'ytick.major.size': 4,
        'ytick.major.pad': 4,
        'ytick.major.width': 1.0
    }
    plt.rcParams.update(params)

    # initialize the report figure and set the title
    fig = plt.figure(layout='tight')
    fig.suptitle('{}-{}-{}, Deployment {:02d}'.format(site.upper(), node.upper(), sensor.upper(), deployment))
    gs = fig.add_gridspec(9, 2)

    # create the plots for each parameter
    ax1 = fig.add_subplot(gs[0, :])
    labels = ['T', 'RH', 'RI']
    colors = ['RoyalBlue', 'DarkOrange', 'ForestGreen']
    yvalues = [1.3, 1, 0.7]
    for i in range(len(availability)):
        if len(availability[i]) > 0:
            ax1.plot(availability[i], availability[i].astype(int) * 0 + yvalues[i], 'o', label=labels[i], color=colors[i])
    ax1.set_xlim(dates)
    ax1.set_ylim([0.5, 1.50])
    ax1.set_yticks(yvalues)
    ax1.set_yticklabels(labels)

    ax2 = fig.add_subplot(gs[1:3, :])
    ax2.plot(merged.time, merged.sea_water_practical_salinity, label='Active', color='RoyalBlue')
    ax2.set_xlim(dates)
    ax2.set_ylim([20, 35])
    ax2.set_yticks([20, 25, 30, 35])
    ax2.set_ylabel('Practical Salinity (psu)')

    ax3 = fig.add_subplot(gs[3:5, 0])
    ax3.plot(merged.time, merged.sea_water_practical_salinity, label='Active', color='RoyalBlue')
    if pre_merged is not None:
        ax3.plot(pre_merged.time, pre_merged.sea_water_practical_salinity, label='Pre', color='grey')
    ax3.scatter(samples['Start Time [UTC]'], samples['CTD Salinity 1 [psu]'], marker='*', label='PSU1')
    ax3.scatter(samples['Start Time [UTC]'], samples['CTD Salinity 2 [psu]'], marker='*', label='PSU2')
    ax3.scatter(samples['Start Time [UTC]'], samples['Discrete Salinity [psu]'], marker='*', label='Discrete')
    zoom = [dates[0] - pd.Timedelta('7D'), dates[0] + pd.Timedelta('7D')]
    ax3.set_xlim(zoom)
    x_fmt = mdates.DateFormatter('%b-%d')
    ax3.xaxis.set_major_formatter(x_fmt)
    ax3.set_xlabel('')
    ax3.set_ylim([20, 35])
    ax3.set_yticks([20, 25, 30, 35])
    ax3.set_ylabel('Practical Salinity (psu)')

    ax4 = fig.add_subplot(gs[3:5, 1])
    ax4.plot(merged.time, merged.sea_water_practical_salinity, label='Active', color='RoyalBlue')
    if post_merged is not None:
        ax4.plot(post_merged.time, post_merged.sea_water_practical_salinity, label='Post', color='grey')
    ax4.scatter(samples['Start Time [UTC]'], samples['CTD Salinity 1 [psu]'], marker='*', label='PSU1')
    ax4.scatter(samples['Start Time [UTC]'], samples['CTD Salinity 2 [psu]'], marker='*', label='PSU2')
    ax4.scatter(samples['Start Time [UTC]'], samples['Discrete Salinity [psu]'], marker='*', label='Discrete')
    zoom = [dates[1] - pd.Timedelta('7D'), dates[1] + pd.Timedelta('7D')]
    ax4.set_xlim(zoom)
    ax4.xaxis.set_major_formatter(x_fmt)
    ax3.set_xlabel('')
    ax4.set_ylim([20, 35])
    ax4.set_yticks([20, 25, 30, 35])

    ax5 = fig.add_subplot(gs[5, :])
    ax5.plot(merged.time, merged.sea_water_electrical_conductivity, label='Active', color='RoyalBlue')
    ax5.set_xlim(dates)
    ax5.set_ylim([2, 5])
    ax5.set_ylabel('Conductivity (S/m)')

    ax6 = fig.add_subplot(gs[6, :])
    ax6.plot(merged.time, merged.sea_water_temperature, label='Active', color='RoyalBlue')
    ax6.set_xlim(dates)
    ax6.set_ylim([5, 20])
    ax6.set_ylabel('Temperature (C)')

    ax7 = fig.add_subplot(gs[7, :])
    ax7.plot(merged.time, merged.sea_water_pressure, label='Active', color='RoyalBlue')
    ax7.set_xlim(dates)
    ax7.set_ylim([26, 32])
    ax7.set_ylabel('Pressure (dbar)')

    ax8 = fig.add_subplot(gs[8, :])
    auto = qartod[0]
    self = qartod[1]
    anno = qartod[2]
    if len(auto) > 0:
        ax8.plot(auto[0].time, auto[0].values, '.', label='S/m')
        ax8.plot(auto[1].time, auto[1].values, '.', label='^oC')
        ax8.plot(auto[2].time, auto[2].values, '.', label='dbar')
        ax8.plot(auto[3].time, auto[3].values, '.', label='psu')
    if len(self) > 0:
        pass  # TODO: add the self-test/older QC-test results to the plot
    if len(anno) > 0:
        pass  # TODO: add the annotation results to the plot

    ax8.set_xlim(dates)
    ax8.set_ylim([0, 5])
    ax8.set_yticks([1, 2, 3, 4])

    # return the figure object
    return fig


def generate_report(site, node, sensor, deployment, cruise_name):
    """
    Load the CTD data for a defined reference designator (using the site, node
    and sensor names to construct the reference designator) collected via the
    telemetered, recovered host and instrument methods and use the data to
    generate a report for HITL QC review.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :param deployment: deployment number to download the data for
    :param cruise_name: name of the deployment cruise(s) to search for
    """
    # load the data from the three delivery methods
    data = combine_delivery_methods(site, node, sensor, deployment)

    # also grab the data from the prior and post deployments
    pre_data = combine_delivery_methods(site, node, sensor, deployment - 1)
    post_data = combine_delivery_methods(site, node, sensor, deployment + 1)

    # create an availability data array for each data stream
    availability = [d.time for d in data if d is not None]

    # get the start and end dates for the deployment
    start, end = get_deployment_dates(site, node, sensor, deployment)
    start = parser.parse(start, ignoretz=True)
    end = parser.parse(end, ignoretz=True)
    dates = [start, end]

    # load the discrete sample data for the deployment cruise
    samples = get_discrete_samples('Endurance', cruise=cruise_name)

    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

    # Convert the flags to QARTOD flags
    codes = {
        None: 0,
        'pass': 1,
        'not_evaluated': 2,
        'suspect': 3,
        'fail': 4,
        'not_operational': 9,
        'not_available': 9,
        'pending_ingest': 0
    }
    annotations['qcFlag'] = annotations['qcFlag'].map(codes).astype('category')

    # limit the annotations to the deployment dates making sure to catch any that might span the deployment, be
    # entirely within the deployment, or have start and/or end dates that fall within the deployment dates
    annotations = annotations[((annotations.beginDate <= start.strftime('%Y-%m-%dT%H:%M:%S')) &
                               (annotations.endDate >= end.strftime('%Y-%m-%dT%H:%M:%S'))) |
                              ((annotations.beginDate >= start.strftime('%Y-%m-%dT%H:%M:%S')) &
                               (annotations.endDate <= end.strftime('%Y-%m-%dT%H:%M:%S'))) |
                              ((annotations.beginDate <= start.strftime('%Y-%m-%dT%H:%M:%S')) &
                               (annotations.endDate >= start.strftime('%Y-%m-%dT%H:%M:%S')) &
                               (annotations.endDate <= end.strftime('%Y-%m-%dT%H:%M:%S'))) |
                              ((annotations.beginDate >= start.strftime('%Y-%m-%dT%H:%M:%S')) &
                               (annotations.beginDate <= end.strftime('%Y-%m-%dT%H:%M:%S')) &
                               (annotations.endDate >= end.strftime('%Y-%m-%dT%H:%M:%S')))]

    # sort the annotations by the beginDate
    annotations = annotations.sort_values(by='beginDate')

    # working through the 3 data streams, reformat the data, add the annotations and apply the QC tests
    for i in range(len(data)):
        if i < 2:
            if data[i]:
                data[i] = ctdbp_datalogger(data[i])
                data[i] = apply_qc_results(data[i], annotations)
            if pre_data[i]:
                pre_data[i] = ctdbp_datalogger(pre_data[i])
                pre_data[i] = apply_qc_results(pre_data[i], annotations)
            if post_data[i]:
                post_data[i] = ctdbp_datalogger(post_data[i])
                post_data[i] = apply_qc_results(post_data[i], annotations)
        else:
            if data[i]:
                data[i] = ctdbp_instrument(data[i])
                data[i] = apply_qc_results(data[i], annotations)
            if pre_data[i]:
                pre_data[i] = ctdbp_instrument(pre_data[i])
                pre_data[i] = apply_qc_results(pre_data[i], annotations)
            if post_data[i]:
                post_data[i] = ctdbp_instrument(post_data[i])
                post_data[i] = apply_qc_results(post_data[i], annotations)

    # combine the data from the three delivery methods into a single dataset
    merged = combine_datasets(data[0], data[1], data[2], None)

    # pull out the QC results for the 4 parameters of interest
    qartod = [
        [merged[x] for x in merged.variables if 'qartod_results' in x],
        [merged[x] for x in merged.variables if 'qc_summary_results' in x],
        [merged[x] for x in merged.variables if 'annotations_qc_results' in x]
    ]

    # now resample the data to hourly intervals (adding a 30-minute offset to the time to center the data)
    merged['time'] = merged['time'] + pd.Timedelta('30Min')
    merged = merged.resample(time='1h', skipna=True).median(dim='time', keep_attrs=True)

    # combine the pre- and post- data from the three delivery methods into a single dataset (pre- and post- data are
    # really used just for plotting purposes)
    pre_merged = combine_datasets(pre_data[0], pre_data[1], pre_data[2], 60)
    post_merged = combine_datasets(post_data[0], post_data[1], post_data[2], 60)

    # trim the pre- and post- data to the deployment dates plus a buffer
    if pre_merged is not None:
        pre_merged = pre_merged.sel(time=slice(start - pd.Timedelta('7D'), start + pd.Timedelta('7D')))
    if post_merged is not None:
        post_merged = post_merged.sel(time=slice(end - pd.Timedelta('7D'), end + pd.Timedelta('7D')))

    # determine the distance to the sampling location for each discrete sample relative to the mooring location
    distance = distance_to_cast(samples, merged.attrs['lat'], merged.attrs['lon'])

    # use the distance to the cast to determine the closest samples to the mooring location within a defined distance
    if site in ['CE01ISSM', 'CE06ISSM']:
        dmin = 0.5  # 0.5 km for the shallow sites
    if site in ['CE02SHSM', 'CE07SHSM']:
        dmin = 1.0  # 1.0 km for the shelf sites
    if site in ['CE040SSM', 'CE09OSSM']:
        dmin = 2.0  # 2.0 km for the offshore sites
    idx = distance[distance < dmin].index
    samples = samples.iloc[idx, :]

    # now find the samples collected within 3 meters of the instrument depth
    if node == 'MFD37':
        drange = 7.5
    else:
        drange = 3.0
    idx = list(samples[(samples['CTD Depth [m]'] - merged['depth'].mean().values).abs() <= drange].index)
    samples = samples.loc[idx]

    # generate the report plot
    fig = plotting(site, node, sensor, deployment, dates, availability, merged, qartod, pre_merged, post_merged, samples)

    # return the figure object and the annotations
    return fig, annotations


def main(argv=None):
    """
    Download the CTDBP data from the Gold Copy THREDDS server and create a
    plot page and annotations table for subsequent HITL QC reviews.
    """
    # set up the input arguments
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    deployment = args.deployment
    cruise_name = args.cruise

    # save the downloaded annotations and plots to a local directory
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/reports/ctdbp')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    fig, annotations = generate_report(site, node, sensor, deployment, cruise_name)

    # save the annotations to a csv file for further processing
    anno_csv = '{}-{}-{}.Deploy{:02d}.quality_annotations.csv'.format(site, node, sensor, deployment)
    annotations.to_csv(os.path.join(out_path, anno_csv), index=False, columns=ANNO_HEADER)

    # save the plot to a pdf file for review
    pdf = '{}-{}-{}.Deploy{:02d}.review_plots.pdf'.format(site, node, sensor, deployment)
    with PdfPages(os.path.join(out_path, pdf)) as pdf_pages:
        pdf_pages.savefig(fig, bbox_inches='tight')

    plt.close(fig)


if __name__ == '__main__':
    main()
