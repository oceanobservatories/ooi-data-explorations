#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Used to download discrete sample data from the OOI Alfresco WebDAV server
"""
import io
import numpy as np
import pandas as pd
import re
import requests

from bs4 import BeautifulSoup


def get_sites(url):
    """
    Get the list of arrays from the OOI Alfresco WebDAV server. The list is
    returned as a list of BeautifulSoup tags that can be further parsed to
    construct URLs to search for the different cruises per array that will
    have discrete sample data.

    :param url: OOI Alfresco WebDAV server URL
    :return sites: list of BeautifulSoup tags containing the array names
    """
    r = requests.get(url, auth=('guest', 'guest'))
    if r.status_code == requests.codes.ok:
        soup = BeautifulSoup(r.text, 'html.parser')
        sites = soup.find_all('a', {'href': re.compile(r'^/alfresco/webdav/OOI/Cabled|Coastal|Global')})
        return sites
    else:
        return None


def get_cruises(url, sites, site_name, cruise=None):
    """
    Get the list of cruises for a specific array from the OOI Alfresco WebDAV
    server. The list is returned as a list of BeautifulSoup tags that can be
    further parsed to construct URLs to search for the different discrete
    sample data files per cruise.

    :param url: OOI Alfresco WebDAV server URL
    :param sites: list of BeautifulSoup tags containing the array names
    :param site_name: name of the array to search for
    :param cruise: optional, name of the cruise(s) to search for (can be
        either a single name or a list of names)
    :return cruises: list of BeautifulSoup tags containing the cruise name(s)
    """
    # if site_name is an empty string, return None (no site to search for)
    if not site_name:
        return None

    # find the site name in the list of sites
    site = [site for site in sites if site_name in site.text]
    if site is None:
        return None

    # create the URL for the site and search for the list of cruises
    a = re.findall(r'<a href=\"(.+?)\">', str(site))[0]
    site_url = url + a + '/Cruise%20Data/'
    r = requests.get(site_url, auth=('guest', 'guest'))
    if r.status_code == requests.codes.ok:
        soup = BeautifulSoup(r.text, 'html.parser')
        if cruise:
            if isinstance(cruise, list):
                cruises = [soup.find_all('a', {'href': re.compile('/Cruise%20Data/.*?' + c + '.+?')})[0] for c in cruise]
            else:
                cruises = soup.find_all('a', {'href': re.compile('/Cruise%20Data/.*?' + cruise + '.+?')})
        else:
            cruises = soup.find_all('a', {'href': re.compile(a + '/Cruise%20Data/')})

        return cruises
    else:
        return None


def get_discrete_samples(site_name, cruise=None):
    """
    Get the discrete sample data for a specific array and cruise from the OOI
    Alfresco WebDAV server. The data is returned as a pandas DataFrame.

    :param site_name: name of the array to search for
    :param cruise: optional, name of the cruise(s) to search for (can be
        either a single name or a list of names)
    :return samples: pandas DataFrame containing the discrete sample data
    """
    # base URL for the OOI Alfresco server
    base_url = 'https://alfresco.oceanobservatories.org'

    # create a listing of the array names in the WebDAV directory
    sites = get_sites(base_url + '/alfresco/webdav/OOI/')

    # use the site name to find the URLs for the different cruises at that site
    cruises = get_cruises(base_url, sites, site_name, cruise)
    if cruises is None:
        # could not find the specified site or cruise(s)
        return None

    # use the list of cruises to find the URLs for the discrete sample CSVs for each cruise and
    # download the data into a pandas DataFrame
    samples = []
    for sample in cruises:
        # create the URL for the cruise discrete sample data directory (some cruises have the data
        # in a slightly different directory)
        sample_url1 = base_url + re.findall(r'<a href=\"(.+?)\">', str(sample))[0] + '/Ship_Data/Water%20Sampling/'
        sample_url2 = base_url + re.findall(r'<a href=\"(.+?)\">', str(sample))[0] + '/Ship%20Data/Water%20Sampling/'

        # search for the discrete sample data files using the two different URLs (note, not all cruises
        # have discrete sample data files)
        r1 = requests.get(sample_url1, auth=('guest', 'guest'))
        r2 = requests.get(sample_url2, auth=('guest', 'guest'))
        if r1.status_code == requests.codes.ok:
            r = r1
        elif r2.status_code == requests.codes.ok:
            r = r2
        else:
            continue

        # find the CSV file for the discrete sample data and create the URL to download the file
        soup = BeautifulSoup(r.text, 'html.parser')
        csv = soup.find_all('a', {'href': re.compile(r'.+?Discrete_Summary\.csv')})
        if not csv:
            continue
        else:
            csv = csv[0]
        csv_url = base_url + re.findall(r'<a href=\"(.+?)\">', str(csv))[0]

        # download the CSV file and read it into a pandas DataFrame appending it to the list of discrete
        # sample data files from the different cruises, if any
        r = requests.get(csv_url, auth=('guest', 'guest'))
        times = ['Start Time [UTC]', 'CTD Bottle Closure Time [UTC]']
        samples.append(pd.read_csv(io.StringIO(r.content.decode('utf-8')), parse_dates=times, date_format='ISO8601'))

    if not samples:
        # could not find any discrete sample data files for the specified site and cruise(s)
        return None

    # concatenate the list of DataFrames into a single DataFrame
    samples = pd.concat(samples, ignore_index=True)

    # replace the fill values with NaNs for all sample columns (ignoring the flag columns
    # and a few other select columns)
    for column in samples.columns:
        if 'Flag' in column:
            # convert the flag columns to a string type
            samples[column] = samples[column].astype(str)
        else:
            if column not in ['Cruise', 'Station', 'Target Asset', 'Cast', 'CTD File',
                              'CTD Bottle Closure Time [UTC]']:
                samples[column] = samples[column].replace(-9999999., value=np.nan)

    return samples
