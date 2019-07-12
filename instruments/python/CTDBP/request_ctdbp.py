#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import yaml

from instruments.python.common import inputs, m2m_collect, m2m_request, dt64_epoch

# load configuration settings
CONFIG = yaml.safe_load(open('instruments\\python\\config.yaml'))


def ctdbp_datalogger(data, burst=False):
    """
    Takes CTDBP data recorded by the data loggers used in the CGSN/EA moorings and cleans up the data set to make
    it more user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that
    are of no use/value. Secondary task is to apply a median average to the burst data collected during some of the
    initial deployments (on some moorings CTDBP was programmed to collect a data point every 10 seconds for 3 minutes).

    :param data: initial ctdbp data set downloaded from OOI via the M2M system
    :param burst: boolean to indicate whether the data was collected in burst mode or not
    :return: cleaned up data set
    """
    # drop some of the variables:
    #   date_time_string == internal_timestamp, redundant so can remove
    #   dcl_controller_timestampe == time, redundant so can remove
    #   id and provenance are not useful in this context. better to obtain directly via M2M request
    # rename some of the variables:
    #   temp is a sloppy name, rename to temperature
    data = data.drop(['date_time_string', 'dcl_controller_timestamp', 'id', 'provenance'])
    #data = data.rename({'temp': 'temperature'})

    # if the data was collected in burst mode, apply a median average to the bursts
    if burst:
        data = data.resample(time='15Min').median()
        data = data.where(~np.isnan(data.temperature), drop=True)

    return data


def ctdbp_instrument(data):
    """
    Takes CTDBP data recorded by internally by the instrument, and cleans up the data set to make it more
    user-friendly. Primary task is renaming the alphabet soup parameter names and dropping some parameters that are
    of no use/value.

    :param data: initial ctdbp data set downloaded from OOI via the M2M system
    :return: cleaned up data set
    """
    return data


def main(argv=None):
    # setup the input arguments
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    method = args.method
    stream = args.stream
    start = args.start
    stop = args.stop

    # request the data
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    ctdbp = m2m_collect(r, '.*CTDBP.*\\.nc')

    # rework the data somewhat to clean-up variables and apply burst averaging
    if method in ['telemetered', 'recovered_host']:
        ctdbp = ctdbp_datalogger(ctdbp)
    else:
        ctdbp = ctdbp_instrument(ctdbp)

    # set the output path to save the downloaded data to
    out_path = CONFIG['base_dir']['m2m_base'] + '\\' + site.lower() + '\\nsif\\ctdbp'
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # create an absolute path and file name combination
    out_file = os.path.abspath(os.path.join(out_path, args.outfile))

    # save the data to disk
    ctdbp.to_netcdf(out_file, mode='w', format='NETCDF4', engine='netcdf4')


if __name__ == '__main__':
    main()
