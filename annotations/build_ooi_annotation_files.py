"""
This script will build a set of files that allow users to sift through OOI annotations. It can be run manually

Requirements:
- OOI credentials must be in a .netrc file located in the home directory.
- SAVE_DIR must be specified, or else the files will be stored in the same directory as the script.
"""

import ast
import csv
from datetime import datetime, timezone
import glob
import json
import logging
import os
import requests
from requests.exceptions import HTTPError
import time
import yaml

SAVE_DIR = os.getcwd() # Declare save directory for files.

# Setup logger for descriptive console messages.
def setup_logger(logger_name: str = 'ooi-annotations',
                 level: int | str = 'DEBUG') -> logging.Logger:
    """
    Set up a logger object for displaying verbose messages to console.

    :param logger_name: The unique logger name to use. Can be shared between modules
    :param level: The logging level to use. Default is 2, which corresponds to DEBUG.
    :return: The configured logging.Logger object.
    """

    logger = logging.getLogger(logger_name)
    logger.propagate = False
    if not logger.handlers:
        logger.setLevel(logging.DEBUG)
        console = logging.StreamHandler()
        console.setLevel(level)

        # Set the logging format.
        dtfmt = '%Y-%m-%dT%H:%M:%S'
        strfmt = f'%(asctime)s.%(msecs)03dZ | %(name)-12s | %(levelname)-8s | %(message)s'
        fmt = logging.Formatter(strfmt, datefmt=dtfmt)
        fmt.converter = time.gmtime

        console.setFormatter(fmt)
        logger.addHandler(console)
    return logger


LOG = setup_logger() # Declare within script scope for use by all functions.


def main():
    LOG.info('Starting build_ooi_annotation_files.py...')

    # Remove old anno files so we don't take up too much space.
    old_anno_files = sorted(glob.glob(os.path.join(SAVE_DIR, "OOI_ANNOTATIONS_UPDATED*.*")))
    if len(old_anno_files) != 0:
        for old_anno_file in old_anno_files:
            os.remove(old_anno_file)
        LOG.info('Removed old annotation files.')

    annos = get_all_annotations()
    to_json(annos, SAVE_DIR)
    to_yaml(annos, SAVE_DIR)
    to_csv(annos, SAVE_DIR)
    LOG.info('build_ooi_annotation_files.py finished.')

    
def get_reference_designators() -> list[str]:
    """
    Obtain a list of OOI reference designators in the form of site - node - instrument.

    :return:
        A list of strings where each string is a unique reference designator.
    """
    LOG.info('Acquiring reference designators...')
    base_url = 'https://ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv'
    base_response = requests.get(base_url)
    rds = []
    if base_response.status_code == requests.codes.ok:
        sites = base_response.json()
        for site in sites:
            site_url = '/'.join([base_url, site]) 
            site_response = requests.get(site_url)
            if site_response.status_code == requests.codes.ok:
                nodes = site_response.json()
                for node in nodes:
                    node_url = '/'.join([site_url, node])
                    node_response = requests.get(node_url)
                    if node_response.status_code == requests.codes.ok:
                        insts = node_response.json()
                        for inst in insts:
                            rd = '-'.join([site, node, inst])
                            rds.append(rd)
                    else:
                        LOG.critical('HTTP Request Error')
                        raise HTTPError(node_response.text)
            else:
                LOG.critical('HTTP Request Error')
                raise HTTPError(site_response.text)
    else:
        LOG.critical('HTTP Request Error')
        raise HTTPError(base_response.text)
    rds = sorted(rds)
    LOG.info('Acquired reference designators.')
    return rds


def get_annotations(refdes: str) -> dict:
    """
    Obtain annotations for a given reference designator.
    This function obtains all annotations created for a given reference designator for all of time.

    :param refdes: 
        A string containing a reference designator in the form of site - node - instrument.
    :return: 
        A dictionary containing the annotation data.
    """
    base_url = 'https://ooinet.oceanobservatories.org/api/m2m/12580/anno/find'
    params = {'beginDT': 0, 
              'endDT': int(datetime(2100,1,1, tzinfo = timezone.utc).timestamp() * 1000),
              'refdes': refdes}
    response = requests.get(base_url, params = params)
    if response.status_code == requests.codes.ok:
        json_data = response.json()
        LOG.info(f'Acquired annotations for {refdes}.')
        return json_data
    else:
        LOG.critical('HTTP Request Error')
        raise HTTPError(response.text)


def int2time(mu_s: int) -> str:
    """
    Convert a returned annotation timestamp (microseconds) to an ISO8601 time string.

    :param us: 
        The microseconds timestamp as an integer.
    :return:
        The timestamp as a string in the format of YYYY-MM-DDTHH:mm:ssZ.
    """
    
    if mu_s is None:
        return None
    elif isinstance(mu_s, int):
        s = mu_s/1000
        dt = datetime.fromtimestamp(s)
        dtstr = dt.strftime('%Y-%m-%dT%H:%M:%SZ')
        return dtstr


def get_parameter_info(parameter_id: int) -> dict:
    """
    Obtain annotations for a given reference designator.
    This function obtains all annotations created for a given reference designator for all of time.

    :param refdes: 
        A string containing a reference designator in the form of site - node - instrument.
    :return: 
        A dictionary containing the annotation data.
    """
    
    base_url = 'https://ooinet.oceanobservatories.org/api/m2m/12575/parameter'
    parameter_url = '/'.join([base_url, str(parameter_id)])
    response = requests.get(parameter_url)
    if response.status_code == requests.codes.ok:
        json_data = response.json()
        return json_data
    else:
        LOG.critical('HTTP Request Error')
        raise HTTPError(response.text)


def param2str(parameters: list[int]) -> str:
    """
    Convert a list of parameter ids to a comma separated string.
    :param parameters: 
        A list of parameter ids.
    :return:
        A stringified list of parameter ids.
    """
    pstr = str(parameters)
    return pstr


def get_param_names(parameters: list[int]) -> str:
    """
    Get parameter names by making secondary requests to the preload API url.

    :param parameters: 
        A list of parameter identifiers.
    :return:
        A stringified list of parameter names.
    """
    if len(parameters) == 0:
        pstr = str(parameters)
    elif len(parameters) == 1:
        param_id = parameters[0]
        param_info = get_parameter_info(param_id)
        if param_info['netcdf_name'] != '':
            pname = param_info['netcdf_name']
        else:
            pname = param_info['name']
        pstr = str([pname])
    else:
        param_infos = [get_parameter_info(param_id) for param_id in parameters]
        names = []
        for param_info in param_infos:
            if param_info['netcdf_name'] != '':
                name = param_info['netcdf_name']
            else:
                name = param_info['name']
            names.append(name) 
        pstr = str(names)
    return pstr


def reformat_annotation(anno: dict) -> dict:
    """
    Reformat a single annotation to be more clear.

    :param anno: 
        An original annotation dictionary as it arrives from the OOI API.
    :return:
        A reformatted dictionary.
    """
    
    anno['site'] = anno.pop('subsite')
    anno['instrument'] = anno.pop('sensor')
    anno['begin_datetime'] = int2time(anno.pop('beginDT'))
    anno['end_datetime'] = int2time(anno.pop('endDT'))
    anno['qc_flag'] = anno.pop('qcFlag')
    anno['exclusion_flag'] = anno.pop('exclusionFlag')
    anno['parameter_ids'] = param2str(anno['parameters']) # For finding unique annotations via list comprehension.
    anno['parameter_names'] = get_param_names(anno.pop('parameters'))
    
    del anno['@class'] # Remove decorator.

    anno = dict(sorted(anno.items())) # Sort keys.
    
    return anno

    
def get_all_annotations() -> list[dict]:
    """
    Get all unique annotations across time and all OOI assets. 
    Output is sorted by begin_datetime and duplicate entries are removed.

    :return:
        A list of dictionaries containing unique annotations.
    """
    LOG.info('Acquiring OOI annotations...')
    rds = get_reference_designators()
    annos = []
    for rd in rds:
        rd_annos = get_annotations(rd)
        if len(rd_annos) == 0:
            continue
        rd_annos = [reformat_annotation(anno) for anno in rd_annos]
        annos.extend(rd_annos)

    # Drop duplicate entries.
    uannos = [dict(s) for s in set(frozenset(d.items()) for d in annos)]

    # Reset dtypes for lists.
    uannos = [{**item, 'parameter_ids': ast.literal_eval(item['parameter_ids'])} for item in uannos]
    uannos = [{**item, 'parameter_names': ast.literal_eval(item['parameter_names'])} for item in uannos]

    # Sort each dictionary alphabetically and then the total list by time.
    uannos = [dict(sorted(d.items())) for d in uannos]
    uannos = sorted(uannos, key = lambda x: x['begin_datetime']) # Sort by time.
    
    LOG.info('All OOI annotations acquired.')    
    
    return uannos


def to_json(_dict: dict, save_directory: os.PathLike) -> None:
    """
    Store a list of annotation dictionaries as JSON.

    :param _dict:
        A list of dictionaries containing annotation data.
    :param save_directory:
        The directory or folder you want to store the annotations file.
    :return: 
        None
    """
    dtstr = datetime.now(timezone.utc).strftime('%Y-%m-%d')
    filename = f"OOI_ANNOTATIONS_UPDATED_{dtstr}.json"
    save_filepath = os.path.join(save_directory, filename)

    with open(save_filepath, 'w') as _file:
        json.dump(_dict, _file, indent = 4)

    if os.path.isfile(save_filepath):
        LOG.info(f"Saved annotations to: {save_filepath}")
    else:
        raise FileNotFoundError(save_filepath)


def to_yaml(_dict: dict, save_directory: os.PathLike) -> None:
    """
    Store a list of annotation dictionaries as YAML.

    :param _dict:
        A list of dictionaries containing annotation data.
    :param save_directory:
        The directory or folder you want to store the annotations file.
    :return: 
        None
    """
    dtstr = datetime.now(timezone.utc).strftime('%Y-%m-%d')
    filename = f"OOI_ANNOTATIONS_UPDATED_{dtstr}.yaml"
    save_filepath = os.path.join(save_directory, filename)

    with open(save_filepath, 'w') as _file:
        yaml.safe_dump(_dict, _file)
        
    if os.path.isfile(save_filepath):
        LOG.info(f"Saved annotations to: {save_filepath}")
    else:
        raise FileNotFoundError(save_filepath)


def to_csv(_dict: dict, save_directory: os.PathLike) -> None:
    """
    Store a list of annotation dictionaries as a CSV.

    :param _dict:
        A list of dictionaries containing annotation data.
    :param save_directory:
        The directory or folder you want to store the annotations file.
    :return: 
        None
    """
        
    dtstr = datetime.now(timezone.utc).strftime('%Y-%m-%d')
    filename = f"OOI_ANNOTATIONS_UPDATED_{dtstr}.csv"
    save_filepath = os.path.join(save_directory, filename)
    
    with open(save_filepath, 'w') as _file:
        fields = _dict[0].keys()
        writer = csv.DictWriter(_file, fieldnames = fields)
        writer.writeheader()
        writer.writerows(_dict)

    if os.path.isfile(save_filepath):
        LOG.info(f"Saved annotations to: {save_filepath}")
    else:
        raise FileNotFoundError(save_filepath)


# Run it.
if __name__ == "__main__":
    main()