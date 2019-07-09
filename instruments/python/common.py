#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@package ooi_dataset_availability.common
@file ooi_dataset_availability/common.py
@author Christopher Wingard
@brief Provides common base classes, definitions and other utilities
"""
import datetime
import glob
import json
import netrc
import os
import re
import requests
import time

import numpy as np
import pandas as pd
import xarray as xr

from bs4 import BeautifulSoup

# setup constants used to access the data from the different M2M interfaces
nc = netrc.netrc()
AUTH = nc.authenticators('ooinet.oceanobservatories.org')
BASE_URL = 'https://ooinet.oceanobservatories.org/api/m2m/'  # base M2M URL
ANNO_URL = '12580/anno/'                                     # Annotations Information
ASSET_URL = '12587/asset/'                                   # Asset and Calibration Information
DEPLOY_URL = '12587/events/deployment/inv/'                  # Deployment Information
SENSOR_URL = '12576/sensor/inv/'                             # Sensor Information
VOCAB_URL = '12586/vocab/inv/'                               # Vocabulary Information


def json2df(infile):
    """
    Read in a JSON formatted data file, and return the results as a panda data frame.
    """
    if not os.path.isfile(infile):
        # if not, return an empty data frame
        print('JSON data file {0} was not found, returning empty data frame'.format(infile))
        return pd.DataFrame()
    else:
        # otherwise, read in the data file
        with open(infile) as jf:
            df = pd.DataFrame(json.load(jf))

        # some of the data files are empty, exit early if so.
        if df.empty:
            print('JSON data file {0} was empty, returning empty data frame'.format(infile))
            return df

        # setup time and the index
        system_epoch = datetime.date(*time.gmtime(0)[0:3])
        ntp_epoch = datetime.date(1900, 1, 1)
        ntp_delta = (system_epoch - ntp_epoch).days * 24 * 3600
        df['time'] = pd.to_datetime(df.port_timestamp - ntp_delta, unit='s')
        df.index = df['time']

        # convert all long integers (int64) to ones acceptable for further processing
        for col in df.columns:
            if df[col].dtype == np.int64:
                df[col] = df[col].astype(np.int32)

        # return the data frame
        return df


def casedpath_unc(path):
    """
    From my friends at Stack Overflow, corrects inconsistencies in the case of path elements. See:

        https://stackoverflow.com/a/35229734

    :param path: best guess at path to glob on while looking for correctly cased path
    :return: corrected path
    """
    unc, p = os.path.splitunc(path)
    r = glob.glob(unc + re.sub(r'([^:/\\])(?=[/\\]|$)', r'[\1]', p))
    return r and r[0] or path


def list_files(url, tag=''):
    """
    Function to create a list of the NetCDF data files in the THREDDS catalog created by a request to the M2M system.

    :param url: URL to user's THREDDS catalog specific to data a data request
    :param tag: regex pattern used to distinguish files of interest
    :return: list of files in the catalog with the full URL path
    """
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    pattern = re.compile(tag)
    return [node.get('href') for node in soup.find_all('a', text=pattern)]


def process_file(catalog_file):
    """
    Function to download one of the NetCDF files, convert to an xarray data set and convert to time as
    the appropriate dimension instead of obs. Drop the extraneous timestamp variables as these were originally not
    intended to be exposed to the user and lead to confusion as to their meaning.

    :param catalog_file: Unique file, referenced by a URL relative to the catalog, to download and convert the data
        file to an xarray data set.
    :return: downloaded data in an xarray dataset.
    """
    dods_url = 'https://opendap.oceanobservatories.org/thredds/dodsC/'
    url = re.sub('catalog.html\\?dataset=', dods_url, catalog_file)
    ds = xr.open_dataset(url).load()  # download the data first, before converting
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.drop(['obs', 'driver_timestamp', 'ingestion_timestamp', 'port_timestamp', 'preferred_timestamp'])

    return ds


def m2m_request(site, node, sensor, method, stream, start=None, stop=None):
    """
    request the data from OOINet via the M2M system

    :param :
    :param :
    :param :
    :param :
    :param :
    :param :
    :param :
    :return:
    """
    # setup the beginning and ending date/time
    if start:
        begin_date = '?beginDT=' + start
    else:
        begin_date = '?beginDT=0'   # use default of 0 to request all data from the beginning of the record

    if stop:
        end_date = '&endDT=' + stop
    else:
        end_date = ''

    options = begin_date + end_date + '&format=application/netcdf'
    r = requests.get(BASE_URL + SENSOR_URL + site + node + sensor + method + stream + options, auth=(AUTH[0], AUTH[2]))
    data = r.json()

    # wait until the request is completed
    check_complete = data['allURLs'][1] + '/status.txt'  # When SOA is actually not that efficient...
    for i in range(1200):
        r = requests.get(check_complete)
        if r.status_code == requests.codes.ok:
            print('request completed')
            break
        else:
            time.sleep(.5)

    # Create a list of the files from the request above using a simple regex as a tag to discriminate the files
    tag = '.*' + sensor[3:8] + '.*\\.nc$'
    files = list_files(data['allURLs'][0], tag)

    # Process the data files found above and concatenate into a single data set
    frames = [process_file(f) for f in files]
    m2m = xr.concat(frames, 'time')
    return m2m
