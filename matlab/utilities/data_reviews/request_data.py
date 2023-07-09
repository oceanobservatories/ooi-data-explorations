#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@package request_data
@file request_data.py
@author Christopher Wingard
@brief Requests, downloads and saves OOI data as a NetCDF file
"""
import getopt
import os
import sys

from  ooi_data_explorations.common import m2m_request, m2m_collect, dt64_epoch

def request_data(site, node, sensor, method, stream, start, stop):
    """
    Request and download data from the OOI M2M system

    :param site:
    :param node:
    :param sensor:
    :param method:
    :param stream:
    :param start:
    :param stop:
    :return data:
    """
    # Request the data (this may take some time).
    r = m2m_request(site, node, sensor, method, stream, start, stop)

    # Use a regex tag to download the sensor data from the THREDDS catalog
    # created by our request.
    #tag = ('.*{}.*\\.nc$'.format(sensor[3:8]))
    tag = ('.*{}.*\\.nc$'.format(stream))
    data = m2m_collect(r, tag)
    return data

def main(argv):
    try:
        opts, args = getopt.getopt(argv, 's:n:r:m:t:b:e:f:', 
            ['site', 'node', 'sensor', 'method', 'stream', 
             'start', 'stop', 'filename'])
    except getopt.GetoptError:
        print('request_data.py -s <site> -n <node> -r <sensor> -m <method> -t <stream> -b <start date> -e <end date> -f <filename>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-s", "--site"):
            site = arg
        if opt in ("-n", "--node"):
            node = arg
        if opt in ("-r", "--sensor"):
            sensor = arg
        if opt in ("-m", "--method"):
            method = arg
        if opt in ("-t", "--stream"):
            stream = arg
        if opt in ("-b", "--start"):
            start = arg
        if opt in ("-e", "--stop"):
            stop = arg
        if opt in ("-f", "--filename"):
            filename = arg

    # request the data
    data = request_data(site, node, sensor, method, stream, start, stop)
    
    # sort by the deployment number and then time
    data = data.sortby(['deployment', 'time'])

    # reset the time record to seconds since 1970 and update the attributes
    data['time'] = dt64_epoch(data.time)
    data.time.attrs['long_name'] = 'Time'
    data.time.attrs['standard_name'] = 'time'
    data.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00 0:00'
    data.time.attrs['calendar'] = 'gregorian'    

    # save the data to disk
    nc_out = os.path.abspath(filename)
    data.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='netcdf4')


if __name__ == "__main__":
   main(sys.argv[1:])
