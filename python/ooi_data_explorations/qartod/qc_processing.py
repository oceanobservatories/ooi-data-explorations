#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Calculates QARTOD test ranges and creates the resulting tables used by
    the OOI QC lookup functions that are used in implementing the QARTOD
    testing.
"""
import numpy as np
import pandas as pd

from ooi_data_explorations.qartod.climatology import Climatology
from ooi_data_explorations.qartod.gross_range import GrossRange


def format_climatology(param, clm, sensor_range, site, node, sensor, stream, source):
    """
    Creates a dictionary object that can later be saved to a CSV formatted
    file for use in the Climatology lookup tables.

    :param param: parameter name of the variable for the calculated user range
    :param clm:
    :param sensor_range:
    :param site: Site designator, extracted from the first part of the reference
        designator
    :param node: Node designator, extracted from the second part of the reference
        designator
    :param sensor: Sensor designator, extracted from the third and fourth part of
        the reference designator
    :param stream: Stream name that contains the data of interest
    :param source: Notes or comments about how the Gross Range values were
        obtained
    :return qc_dict: dictionary with the sensor and user gross range values
        added in the formatting expected by the QC lookup
    """
    # create the lookup dictionary
    var_explained = clm.regression['variance_explained']
    qc_dict = {
        'subsite': site,
        'node': node,
        'sensor': sensor,
        'stream': stream,
        'parameters': {'inp': param, 'tinp': 'time', 'zinp': 'None'},
        'climatologyTable': 'climatology_tables/{}-{}-{}-{}.csv'.format(site, node, sensor, param),
        'source': source + '. The variance explained by the climatological model is {:.1%}.'.format(var_explained[0])
    }

    # create the climatology table
    header_str = ''
    value_str = '[0, 0]'
    for idx, mu in enumerate(clm.monthly_fit):
        # use the index number to create the header row
        header_str += ',[{}, {}]'.format(idx+1, idx+1)

        # calculate the climatological ranges
        cmin = mu - clm.monthly_std.values[idx] * 3
        if cmin < sensor_range[0]:
            cmin = sensor_range[0]

        cmax = mu + clm.monthly_std.values[idx] * 3
        if cmax > sensor_range[1]:
            cmax = sensor_range[1]

        # append the data to ranges
        value_str += ',[{:.2f}, {:.2f}]'.format(cmin, cmax)

    clm_table = header_str + '\n' + value_str

    return qc_dict, clm_table


def process_climatology(ds, params, sensor_range, **kwargs):
    """

    :param ds:
    :param params:
    :param sensor_range:
    :param kwargs:
    :return:
    """
    # process the optional keyword arguments
    site = kwargs.get('site')
    node = kwargs.get('node')
    sensor = kwargs.get('sensor')
    stream = kwargs.get('stream')
    source = kwargs.get('source')

    # initialize the Climatology class
    clm = Climatology()

    # create an empty panda dataframe and list to hold the results
    clm_lookup = pd.DataFrame()
    clm_tables = []

    # loop through the parameter(s) of interest
    sensor_range = np.atleast_2d(sensor_range).tolist()
    for idx, param in enumerate(params):
        # calculate the 2-cycle climatology for the parameter of interest
        clm.fit(ds, param)

        # create the formatted dictionary for the lookup tables
        qc_dict, clm_table = format_climatology(param, clm, sensor_range[idx], site, node, sensor, stream, source)

        # append the dictionary to the dataframe and the table to the list
        clm_lookup = clm_lookup.append(qc_dict, ignore_index=True)
        clm_tables.append(clm_table)

    # return the results
    return clm_lookup, clm_tables


def format_gross_range(param, sensor_range, user_range, site, node, sensor, stream, source):
    """
    Creates a dictionary object that can later be saved to a CSV formatted
    file for use in the Gross Range lookup tables.

    :param param: parameter name of the variable for the calculated user range
    :param sensor_range: default sensor, or fail range, usually referenced
        from the vendor documentation
    :param user_range: user range, or sensor range, calculated from the data
    :param site: Site designator, extracted from the first part of the reference
        designator
    :param node: Node designator, extracted from the second part of the reference
        designator
    :param sensor: Sensor designator, extracted from the third and fourth part of
        the reference designator
    :param stream: Stream name that contains the data of interest
    :param source: Notes or comments about how the Gross Range values were
        obtained
    :return qc_dict: dictionary with the sensor and user gross range values
        added in the formatting expected by the QC lookup tables
    """
    # create the dictionary
    qc_dict = {
        'subsite': site,
        'node': node,
        'sensor': sensor,
        'stream': stream,
        'parameter': {
            'inp': param
        },
        'qcConfig': {
             'qartod': {
                 'gross_range_test': {
                     'suspect_span': ['{}'.format(user_range[0]), '{}'.format(user_range[1])],
                     'fail_span': ['{}'.format(sensor_range[0]), '{}'.format(sensor_range[1])]
                 }
             }
         },
        'source': source
    }
    return qc_dict


def process_gross_range(ds, params, sensor_range, **kwargs):
    """

    :param ds:
    :param params:
    :param sensor_range:
    :param kwargs:
    :return:
    """
    # process the optional keyword arguments
    site = kwargs.get('site')
    node = kwargs.get('node')
    sensor = kwargs.get('sensor')
    stream = kwargs.get('stream')
    source = kwargs.get('source')

    # create an empty pandas dataframe to hold the results
    gross_range = pd.DataFrame()

    # loop through the parameter(s) of interest
    sensor_range = np.atleast_2d(sensor_range).tolist()
    for idx, param in enumerate(params):
        # calculate the user range
        gr = GrossRange(fail_min=sensor_range[idx][0], fail_max=sensor_range[idx][1])
        gr.fit(ds, param, 3)
        user_range = [gr.suspect_min, gr.suspect_max]
        # create the formatted dictionary
        qc_dict = format_gross_range(param, sensor_range[idx], user_range, site, node, sensor, stream, source)
        # append the dictionary to the dataframe
        gross_range = gross_range.append(qc_dict, ignore_index=True)

    # return the results
    return gross_range
