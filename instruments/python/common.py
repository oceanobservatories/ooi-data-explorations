#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Provides common base classes, definitions and other utilities used to access and download OOI data via the
    M2M interface
"""
import argparse
import netrc
import re
import requests
import sys
import time
import xarray as xr

from bs4 import BeautifulSoup
from tqdm import tqdm

# setup constants used to access the data from the different M2M interfaces
BASE_URL = 'https://ooinet.oceanobservatories.org/api/m2m/'  # base M2M URL
ANNO_URL = '12580/anno/'                                     # Annotation Information
ASSET_URL = '12587/asset/'                                   # Asset and Calibration Information
DEPLOY_URL = '12587/events/deployment/inv/'                  # Deployment Information
SENSOR_URL = '12576/sensor/inv/'                             # Sensor Information
VOCAB_URL = '12586/vocab/inv/'                               # Vocabulary Information

# setup access credentials
nc = netrc.netrc()  # best option is user has their own account
AUTH = nc.authenticators('ooinet.oceanobservatories.org')
if AUTH is None:
    # Use our default user credentials.
    AUTH = ['OOIAPI-853A3LA6QI3L62', 'WYAN89W5X4Z0QZ']
    raise Warning('Using default access credentials for the OOI Coastal Endurance Data Team. ' +
                  'User is encouraged to setup their own access credentials in a netrc file ' +
                  'maintained in their home directory using ooinet.oceanobservatories.org as the ' +
                  'machine name.')


def list_nodes(site):
    """
    Based on the site name, list the nodes that are available.

    :param site: Site name to query
    :return: List of the the available nodes for this site
    """
    r = requests.get(BASE_URL + DEPLOY_URL + site, auth=(AUTH[0], AUTH[2]))
    return r.json()


def list_sensors(site, node):
    """
    Based on the site and node name, list the sensors that are available.

    :param site: Site name to query
    :param node: Node name to query
    :return: list of the the available sensors for this site and node
    """
    r = requests.get(BASE_URL + DEPLOY_URL + site + '/' + node, auth=(AUTH[0], AUTH[2]))
    return r.json()


def list_methods(site, node, sensor):
    """
    Based on the site, node and sensor name, list the data delivery methods that are available.

    :param site: Site name to query
    :param node: Node name to query
    :param sensor: Sensor name to query
    :return: list of the data delivery methods associated with this site, node and sensor
    """
    r = requests.get(BASE_URL + SENSOR_URL + site + '/' + node + '/' + sensor, auth=(AUTH[0], AUTH[2]))
    return r.json()


def list_streams(site, node, sensor, method):
    """
    Based on the site, node and sensor name and the data delivery method, list the data streams that are available.

    :param site: Site name to query
    :param node: Node name to query
    :param sensor: Sensor name to query
    :param method: Data delivery method to query
    :return: list of the data streams associated with this site, node, sensor and data delivery method
    """
    # Determine the streams associated with the delivery method available for this sensor
    r = requests.get(BASE_URL + SENSOR_URL + site + '/' + node + '/' + sensor + '/' + method, auth=(AUTH[0], AUTH[2]))
    return r.json()


def list_deployments(site, node, sensor):
    """
    Based on the site, node and sensor name, list the deployment numbers that are available.

    :param site: Site name to query
    :param node: Node name to query
    :param sensor: Sensor name to query
    :return: list of the deployments of this site, node and sensor combination
    """
    r = requests.get(BASE_URL + DEPLOY_URL + site + '/' + node + '/' + sensor, auth=(AUTH[0], AUTH[2]))
    return r.json()


def deployment_dates(site, node, sensor, deploy):
    """
    Based on the site, node and sensor names and the deployment number, determine the start and end times for a
    deployment.
    :param site:
    :param node:
    :param sensor:
    :param deploy:
    :return:
    """
    # request deployment metadata
    r = requests.get(BASE_URL + DEPLOY_URL + site + '/' + node + '/' + sensor + '/' + str(deploy), auth=(AUTH[0], AUTH[2]))
    data = r.json()

    # use the metadata to extract the start and end times for the deployment
    start = time.strftime('%Y-%m-%dT%H:%M:%S.000Z', time.gmtime(data[0]['eventStartTime'] / 1000.))
    if data[0]['eventStopTime']:
        # check to see if there is a stop time for the deployment, if so use it ...
        stop = time.strftime('%Y-%m-%dT%H:%M:%S.000Z', time.gmtime(data[0]['eventStoptTime'] / 1000.))
    else:
        # ... otherwise use the current time as this is an active deployment
        stop = time.strftime('%Y-%m-%dT%H:%M:%S.000Z', time.gmtime(time.time()))

    return start, stop


def list_files(url, tag=''):
    """
    Function to create a list of the NetCDF data files in the THREDDS catalog created by a request to the M2M system.

    :param url: URL to user's THREDDS catalog specific to a data request
    :param tag: regex pattern used to distinguish files of interest
    :return: list of files in the catalog with the URL path set relative to the catalog
    """
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    pattern = re.compile(tag)
    return [node.get('href') for node in soup.find_all('a', text=pattern)]


def process_file(catalog_file):
    """
    Function to download one of the NetCDF files as an xarray data set, convert to time as the appropriate dimension
    instead of obs, and drop the extraneous timestamp variables (these were originally not intended to be exposed to
    users and lead to confusion as to their meaning).

    :param catalog_file: Unique file, referenced by a URL relative to the catalog, to download and convert the data
        file to an xarray data set.
    :return: downloaded data in an xarray dataset.
    """
    dods_url = 'https://opendap.oceanobservatories.org/thredds/dodsC/'
    url = re.sub('catalog.html\\?dataset=', dods_url, catalog_file)
    ds = xr.open_dataset(url).load()  # fully download the data first, before converting
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.drop(['obs', 'driver_timestamp', 'ingestion_timestamp', 'port_timestamp', 'preferred_timestamp'])

    # convert datetime values from nanoseconds to seconds since 1970-01-01
    ds['time'] = dt64_epoch(ds.time)
    ds['internal_timestamp'] = dt64_epoch(ds.internal_timestamp)

    # resetting cdm_data_type from Point to Station, more accurate for time series datasets
    ds.attrs['cdm_data_type'] = 'Station'

    return ds


def m2m_request(site, node, sensor, method, stream, start=None, stop=None):
    """
    Request data from OOINet for a particular instrument (as defined by the reference designator), delivery method
    and stream name via the M2M system. Optionally, can bound the data with a beginning and ending date and time range.

    :param site: Site designator, extracted from the first part of the reference designator
    :param node: Node designator, extracted from the second part of the reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part of the reference designator
    :param method: Delivery method for the data (either telemetered, recovered_host, or recovered_inst)
    :param stream: Stream name that contains the data of interest
    :param start: Start time for data request (Optional, default is beginning of record)
    :param stop: Stop time for data request (Optional, default is through the end of the record)
    :return: The results of the data request detailing where the data is located for download
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
    r = requests.get(BASE_URL + SENSOR_URL + site + '/' + node + '/' + sensor + '/' + method + '/' + stream + options,
                     auth=(AUTH[0], AUTH[2]))
    data = r.json()

    # wait until the request is completed
    check_complete = data['allURLs'][1] + '/status.txt'
    with tqdm(total=1200, desc='Processing Request') as bar:
        for i in range(1200):
            r = requests.get(check_complete)
            bar.update(1)
            if r.status_code == requests.codes.ok:
                bar.n = 1200
                bar.last_print_n = 1200
                bar.refresh()
                print('\nrequest completed in %f minutes\n' % elapsed)
                break
            else:
                time.sleep(.5)
                elapsed = (i * 0.5) / 60

    return data


def m2m_collect(data, tag=''):
    """
    Use a regex tag combined with the results of the M2M data request to collect the data from the THREDDS catalog.
    Collected data is gathered into an xarray dataset for further processing.

    :param data: JSON object returned from M2M data request with details on where the data is to be found for download
    :param tag: regex tag to use in discriminating the data files, so we only collect the correct ones
    :return: the collected data as an xarray dataset
    """
    # Create a list of the files from the request above using a simple regex as a tag to discriminate the files
    files = list_files(data['allURLs'][0], tag)

    # Process the data files found above and concatenate into a single data set
    frames = [process_file(f) for f in files]
    m2m = xr.concat(frames, 'time')
    return m2m


def dt64_epoch(dt64):
    """
    Convert a panda or xarray date/time value represented as a datetime64 object (nanoseconds since 1970) to a float,
    representing an epoch time stamp (seconds since 1970-01-01).

    :param dt64: panda or xarray datatime64 object
    :return epts: epoch time as seconds since 1970-01-01
    """
    epts = dt64.values.astype(float) / 10.0 ** 9
    return epts


def inputs(argv=None):
    """
    Sets the main input arguments that would be passed by the M2M requesting module.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize argument parser
    parser = argparse.ArgumentParser(description="""Request and obtain data from the OOI M2M system""",
                                     epilog="""Obtain data from the OOI M2M system""")

    # assign input arguments.
    parser.add_argument("-s", "--site", dest="site", type=str, required=True)
    parser.add_argument("-n", "--node", dest="node", type=str, required=True)
    parser.add_argument("-sn", "--sensor", dest="sensor", type=str, required=True)
    parser.add_argument("-mt", "--method", dest="method", type=str, required=True)
    parser.add_argument("-st", "--stream", dest="stream", type=str, required=True)
    parser.add_argument("-bt", "--beginDT", dest="start", type=str, required=False)
    parser.add_argument("-et", "--endDT", dest="stop", type=str, required=False)
    parser.add_argument("-et", "--endDT", dest="stop", type=str, required=False)
    parser.add_argument("-o", "--outfile", dest="outfile", type=str, required=False)

    # parse the input arguments and create a parser object
    args = parser.parse_args(argv)

    return args
