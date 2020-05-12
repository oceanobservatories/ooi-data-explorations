#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Creates the YAML file with the information needed to construct URLs for M2M requests
"""
import requests

from ooi_data_explorations.common import AUTH, BASE_URL
from ooi_data_explorations.common import list_nodes, list_sensors, list_methods, list_streams

# Organize the nodes into assemblies common across all the arrays. Helps to better organize the data, taking all the
# different names used in OOI and collapsing them down to a more coherent, and logical grouping. This list is
# deliberately non-inclusive of all possible nodes, as some are either specific to control systems or have data not
# readily accessed through the M2M system.
ASSEMBLY = {
    'surface_buoy': ['SBD11', 'SBD12', 'SBD17'],
    'midwater': ['PC01A', 'PC01B', 'PC03A', 'RID16', 'RID26', 'RID27', 'RII01', 'RII11', 'RIM01', 'RIS01'],
    'seafloor': [
        'LJ01A', 'LJ01B', 'LJ01C', 'LJ01D', 'LJ03A', 'MFD35', 'MFD37', 'MJ01A', 'MJ01B', 'MJ03A', 'MJ03B',
        'MJ03C', 'MJ03D', 'MJ03E', 'MJ03F', 'PN03B'],
    'glider': [
        'GL247', 'GL276', 'GL311', 'GL312', 'GL319', 'GL320', 'GL326', 'GL327', 'GL335', 'GL336', 'GL339',
        'GL340', 'GL361', 'GL362', 'GL363', 'GL364', 'GL365', 'GL374', 'GL375', 'GL376', 'GL379', 'GL380',
        'GL381', 'GL382', 'GL383', 'GL384', 'GL386', 'GL387', 'GL388', 'GL389', 'GL453', 'GL469', 'GL470',
        'GL477', 'GL478', 'GL484', 'GL485', 'GL486', 'GL493', 'GL494', 'GL495', 'GL496', 'GL514', 'GL523',
        'GL524', 'GL525', 'GL537', 'GL538', 'GL559', 'GL560', 'GL561', 'PG514', 'PG515', 'PG528', 'PG562',
        'PG563', 'PG564', 'PG565', 'PG566', 'PG575', 'PG578', 'PG580', 'PG583'
    ],
    'auv': ['A6263', 'A6264'],
    'profiler': ['WFP01', 'WFP02', 'WFP03', 'SP001', 'SF01A', 'SF01B', 'SF03A', 'DP01A', 'DP01B', 'DP03A']
}

# Assemblies can be further categorized by a sub-assembly type if needed. Otherwise, the code and assembly type are 
# the same. This helps to differentiate times where more than one assembly type is present at a site and we need to 
# further distinguish between them. 
SUBASSEMBLY = {
    'nsif': ['RID16', 'RID26', 'RID27'],
    'riser': ['RII01', 'RII11', 'RIM01', 'RIS01']
}

# List of allowed methods (there are some called bad_*, ignoring those)
METHODS = ['streamed', 'telemetered', 'recovered_host', 'recovered_inst', 'recovered_cspp', 'recovered_wfp']

# Sensor codes to exclude. These are either not science sensors (e.g. 00-DCLENG000), or not really accessible from the
# M2M system (e.g. cameras, hydrophones or seismometers).
SENSOR_EXCLUDES = [
    'AOABP', 'CAM', 'COVIS', 'ENG', 'FLOBN', 'HPIES', 'HYD', 'MASSP', 'OBS', 'OSMOI', 'OVRSR', 'PPSDNA', 'QNTSR',
    'RASFL', 'SCPRA', 'SCTAA', 'ZPLSCB',
]

# Several streams were defined within OOI that are beyond core science needs. The terms below filter out those streams.
STREAM_EXCLUDES = [
    'metadata', 'diagnostics', 'abc_power', 'dcl_power', 'config', 'wave_burst', 'imodem_status',
    'imodem_start_time', 'cspp_dark_instrument', 'instrument_blank', 'adcp_engineering', 'adcp_ancillary',
    'adcp_transmit', 'dcl_fourier', 'dcl_mean_directional', 'dcl_motion', 'dcl_non_directional', 'echogram_data',
    'calibration', 'status', 'settings', 'hardware', 'adcp_pd0_beam_parsed', 'identification_string', 'clock_data',
    'battery_voltage', 'thermistor_voltage', 'dev1_data', 'data_record_cal', 'dark_sample', 'imodem_control',
    'imodem_power', 'offset', 'adcp_velocity_glider', 'nutnr_b_dcl_dark', 'nutnr_b_metadata', 'nutnr_b_dark',
    'prest_event_counter', 'prest_reference_oscillator', 'vel3d_b_engineering', 'tmpsf_engineering',
    'vadcp_5thbeam_pd0_beam_parsed', 'nutnr_a_test'
]

# grab the entire OOI vocabulary dictionary
def get_vocabulary():
    '''
    Based on the site, node and sensor name download the vocabulary record defining this sensor.

    :param site: Site name to query
    :return: json object with the site-node-sensor specific vocabulary
    '''
    r = requests.get(BASE_URL + '12586/vocab', auth=(AUTH[0], AUTH[2]))
    if r.status_code == requests.codes.ok:
        return r.json()
    else:
        return None

VOCAB = get_vocabulary()

def add_site(site, url_file):
    '''
    For an OOI site, assemble a curated list of all the instruments and data streams that are available for data
    explorations. This file will create a YAML file per site. Additional HITL work is required to further clean-up and
    check the list, and pruning streams and methods down to the core set of science sensors. The results of this work
    is used to create the m2m_urls.yml file used in this python toolbox.

    :param site: OOI 8 character site designation (e.g. CE01ISSM)
    :param url_file: data file to save the results in.
    :return: None, creates the data file.
    '''
    sp = '  '  # create an explicit two-space indent string for the YAML file
    # open the YAML file to store the results
    with open(url_file, 'w') as f:
        # create the header portion of the YAML, getting the array and site names from the API
        site = site.upper()
        for site_vocab in VOCAB:
            if site_vocab.get('refdes') == site:
                break

        f.write('---\n{site}:\n'.format(site=site))
        f.write('{sp}array: {array}\n'.format(sp=sp*1, array=site_vocab['tocL1']))
        f.write('{sp}name: {name}\n'.format(sp=sp*1, name=site_vocab['tocL2']))
        f.write('{sp}assembly:\n'.format(sp=sp*1))

        # create a list of nodes for this site
        nodes = list_nodes(site)

        # for each node, if it is one of interest as defined above, create an assembly entry
        for k, v in ASSEMBLY.items():
            # find the nodes that correspond to this assembly
            assembly = k
            node = sorted(set(v) & set(nodes))

            # skip if we don't have any nodes in this assembly
            if not node:
                continue

            # if we have nodes, create the assembly entry
            for node_vocab in VOCAB:
                if node_vocab.get('refdes') == site + '-' + node[0]:
                    break

            # for each node, create the list of sensors
            for n in node:
                f.write('{sp}- type: {assembly}\n'.format(sp=sp * 2, assembly=assembly))
                f.write('{sp}name: {name}\n'.format(sp=sp * 3, name=node_vocab['tocL3']))
                # if we need to further distinguish, add the assembly code name
                for k, v in SUBASSEMBLY.items():
                     if n in v:
                         f.write('{sp}subassembly: {name}\n'.format(sp=sp * 3, name=k))
                f.write('{sp}instrument:\n'.format(sp=sp * 3))
                sensors = list_sensors(site, n)
                sensors = filter_stream(sensors, SENSOR_EXCLUDES)   # remove sensors of no interest

                if not sensors:
                    continue

                for sensor in sensors:
                    for sensor_vocab in VOCAB:
                        if sensor_vocab.get('refdes') == site + '-' + n + '-' + sensor:
                            break

                    instrument = (sensor[3:8]).lower()
                    f.write('{sp}- class: {instrument}\n'.format(sp=sp*4, instrument=instrument))
                    f.write('{sp}instrument_name: {name}\n'.format(sp=sp*5, name=sensor_vocab['instrument']))
                    f.write('{sp}instrument_model: {model}\n'.format(sp=sp*5, model=sensor_vocab['model']))
                    f.write('{sp}instrument_manufacturer: {manu}\n'.format(sp=sp*5, manu=sensor_vocab['manufacturer']))
                    f.write('{sp}mindepth: {mindepth}\n'.format(sp=sp*5, mindepth=sensor_vocab['mindepth']))
                    f.write('{sp}maxdepth: {maxdepth}\n'.format(sp=sp*5, maxdepth=sensor_vocab['maxdepth']))
                    f.write('{sp}node: {node}\n'.format(sp=sp*5, node=n))
                    f.write('{sp}sensor: {sensor}\n'.format(sp=sp*5, sensor=sensor))
                    f.write('{sp}stream:\n'.format(sp=sp*5))

                    methods = list_methods(site, n, sensor)
                    if not methods:
                        f.write('{sp}unknown: null\n'.format(sp=sp*6,))
                        continue

                    for method in methods:
                        if method in METHODS:
                            streams = list_streams(site, n, sensor, method)
                            streams = filter_stream(streams, STREAM_EXCLUDES)
                            if len(streams) == 1:
                                f.write('{sp}{method}: {streams}\n'.format(sp=sp*6, method=method,
                                                                               streams=streams[0]))
                            else:
                                f.write('{sp}{method}:\n'.format(sp=sp*6, method=method))
                                for stream in streams:
                                    f.write('{sp}- {stream}\n'.format(sp=sp*7, stream=stream))

def filter_stream(streams, excludes):
    '''
    Uses a list of keywords to remove sensors or streams from the list returned by OOI Net.

    :param streams: list of sensor or streams returned from OOI Net
    :param excludes: list of keywords to use in pruning the list
    :return: a cleaned, pruned list
    '''
    clean = []
    for stream in streams:
        if not any(sub in stream for sub in excludes):
            clean.append(stream)

    return clean
