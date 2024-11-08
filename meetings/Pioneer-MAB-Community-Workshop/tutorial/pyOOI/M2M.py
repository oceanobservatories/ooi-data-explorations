import os
import re
import sys
import time
import netrc
import requests
import datetime
import numpy as np
import pandas as pd
import xarray as xr
from queue import Queue
from bs4 import BeautifulSoup
from urllib.request import urlretrieve

# Add paths
sys.path.append('../')

# Import shared utilies
from pyOOI.utils import ntp_seconds_to_datetime, convert_time, unix_epoch_time
from pyOOI.Download import setup_download_dir, download_file, DownloadWorker

# Initialize credentials
try:
    nrc = netrc.netrc()
    AUTH = nrc.authenticators('ooinet.oceanobservatories.org')
    login, password = AUTH[0], AUTH[2]
    if AUTH is None:
        raise RuntimeError(
            'No entry found for machine ``ooinet.oceanobservatories.org`` in the .netrc file')
except FileNotFoundError as e:
    raise OSError(e, os.strerror(e.errno), os.path.expanduser('~'))

# Set the URL endpoints for ooinet.oceanobservatories.org
URLS = {
    'data': 'https://ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv',
    'anno': 'https://ooinet.oceanobservatories.org/api/m2m/12580/anno/find',
    'vocab': 'https://ooinet.oceanobservatories.org/api/m2m/12586/vocab/inv',
    'asset': 'https://ooinet.oceanobservatories.org/api/m2m/12587',
    'deploy': 'https://ooinet.oceanobservatories.org/api/m2m/12587/events/deployment/inv',
    'preload': 'https://ooinet.oceanobservatories.org/api/m2m/12575/parameter',
    'cal': 'https://ooinet.oceanobservatories.org/api/m2m/12587/asset/cal',
    'fileServer':"https://opendap.oceanobservatories.org/thredds/fileServer/",
    'dodsC': "https://opendap.oceanobservatories.org/thredds/dodsC/",
    "goldCopy": "https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/",
    "goldCopy_fileServer": "https://thredds.dataexplorer.oceanobservatories.org/thredds/fileServer/",
    "goldCopy_dodsC": "https://thredds.dataexplorer.oceanobservatories.org/thredds/dodsC/"
}

# Set the connections
SESSION = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=0)
SESSION.mount("https://", adapter)


def get_api(url, params=None):
    """Function which gets OOINet API endpoint"""
    r = SESSION.get(url, params=params, auth=(login, password))
    if r.status_code == requests.codes.ok:
        api = r.json()
    else:
        print(r.reason)
        api = None
    return api


def get_datasets(search_url, datasets=pd.DataFrame(), **kwargs):
    """Search OOINet for available datasets for a url."""
    # Check if the method is attached to the url
    flag = ("inv" == search_url.split("/")[-4])
    # inst = re.search("[0-9]{2}-[023A-Z]{6}[0-9]{3}", search_url)
    # inst = re.search("[0-9]{2}-", search_url)

    # This means you are at the end-point
    if flag is True:
        # Get the reference designator info
        array, node, instrument = search_url.split("/")[-3:]
        refdes = "-".join((array, node, instrument))

        # Get the available deployments
        deploy_url = "/".join((URLS["deploy"], array, node,
                               instrument))
        deployments = get_api(deploy_url)

        # Put the data into a dictionary
        info = {
            "array": array,
            "node": node,
            "instrument": instrument,
            "refdes": refdes,
            "url": search_url,
            "deployments": [deployments],
        }
        
        # Convert it to a dataframe
        info_df = pd.DataFrame(info)
        
        # add the dictionary to the dataframe
        datasets = pd.concat([datasets, info_df])

    else:
        endpoints = get_api(search_url)
                   
        while len(endpoints) > 0:

            # Get one endpoint
            new_endpoint = endpoints.pop()

            # Build the new request url
            new_search_url = "/".join((search_url, new_endpoint))

            # Get the datasets for the new given endpoint
            datasets = get_datasets(new_search_url, datasets)

    # Once recursion is done, return the datasets
    return datasets


def search_datasets(array=None, node=None, instrument=None,
                    English_names=False):
    """Search OOINet for datasets.

    Parameters
    ----------
    array: (str)
        OOI abbreviation for a particular buoy on an array (e.g. Pioneer
        Central Surface Mooring = CP01CNSM)
    node: (str)
        Partial or full OOI abbreviation for a node on a buoy to search for
        (e.g. Multi-Function Node = MFD)
    instrument: (str)
        Partial or full OOI abbreviation for a particular instrument type
        to search for (e.g. CTD)
    English_names: (bool)
        Set to True if the descriptive names associated with the given
        array/node/instrument are wanted.

    Returns
    -------
    datasets: (pandas.DataFrame)
        A dataframe of all the OOI datasets which match the given search
        terms. If no search terms are entered, will return every dataset
        available in OOINet (slow).
    """
    # Build the request url
    dataset_url = f'{URLS["data"]}/{array}/{node}/{instrument}'

    # Truncate the url at the first "none"
    dataset_url = dataset_url[:dataset_url.find("None")-1]

    print(f"Searching {dataset_url}")
    # Get the datasets
    datasets = get_datasets(dataset_url)

    # Now, it node is not None, can filter on that
    if node is not None:
        mask = datasets["node"].apply(lambda x: True if node
                                      in x else False)
        datasets = datasets[mask]

    # If instrument is not None
    if instrument is not None:
        mask = datasets["instrument"].apply(lambda x: True if instrument
                                            in x else False)
        datasets = datasets[mask]

    # Check if they want the English names for the associated datasets
    if English_names:
        vocab = {
            "refdes": [],
            "array_name": [],
            "node_name": [],
            "instrument_name": []
        }

        # Iterate through the given reference designators
        for refdes in datasets["refdes"]:
            # Request the vocab for the given reference designator
            refdes_vocab = get_vocab(refdes)

            # Check if it returns an empty dataframe - then fill with NaNs
            if len(refdes_vocab) == 0:
                vocab["refdes"].append(refdes)
                vocab["array_name"].append("No record")
                vocab["node_name"].append("No record")
                vocab["instrument_name"].append("No record")

            # If it isn't empty - Parse the refdes-specific vocab
            else:
                vocab["refdes"].append(refdes)
                vocab["array_name"].append(
                    refdes_vocab["tocL1"].iloc[0] + " " +
                    refdes_vocab["tocL2"].iloc[0])
                vocab["node_name"].append(refdes_vocab["tocL3"].iloc[0])
                vocab["instrument_name"].append(
                    refdes_vocab["instrument"].iloc[0])

        # Merge the results with the datasets
        vocab = pd.DataFrame(vocab)
        datasets = datasets.merge(vocab, left_on="refdes",
                                  right_on="refdes")
        # Sort the datasets
        columns = ["array", "array_name", "node", "node_name",
                   "instrument", "instrument_name", "refdes", "url",
                   "deployments"]
        datasets = datasets[columns]

    return datasets


def get_vocab(refdes):
    """Get OOI vocabulary.

    Return the OOI vocabulary for a given url endpoint. The vocab results
    contains info about the reference designator, names of the

    Parameters
    ----------
    refdes: (str)
        The reference designator for the instrument for which to request
        vocab information.

    Returns
    -------
    results: (pandas.DataFrame)
        A table of the vocab information for the given reference
        designator.
    """
    # First, construct the vocab request url
    array, node, instrument = refdes.split("-", 2)
    vocab_url = "/".join((URLS["vocab"], array, node, instrument))

    # Next, get the vocab data
    data = get_api(vocab_url)

    # Put the returned vocab data into a pandas dataframe
    vocab = pd.DataFrame()
    vocab = pd.concat([vocab, pd.DataFrame(data)])

    # Finally, return the results
    return vocab


def reformat_calInfo(calInfo, deployNum, uid):
    """Internal method to reformat calibration info"""
    calDataFrame = pd.DataFrame()
    # Reformat the calibration info
    for cal in calInfo["calibration"]:
        for calData in cal["calData"]:
            df = pd.DataFrame({
                "deploymentNumber": [int(deployNum)],
                "uid": [uid],
                "calCoef": list(calData["eventName"]),
                "value": list(calData["value"]),
                "calFile": list(calData["dataSource"])
            })
            calDataFrame = pd.concat([calDataFrame, df], ignore_index=True)
    return calDataFrame


def get_calibrations_by_refdes(refdes, deployments):
    """Get calibrations for deployments for a given reference designator.

    Return the effective calibrations for a given reference designator
    for all of the deployments.

    Parameters
    ----------
    refdes: (str)
        The reference designator for the instrument for which to request
        vocab information.
    deployments: (pandas.DataFrame)
        A dataframe with the given deployments for a reference designator
        as returned by the get_deployments method.

    Returns
    -------
    calibrations: (pandas.DataFrame)
        A table of the calibration information for the given reference
        designator for all of the deployments specified by the
        deployments dataframe.
    """
    # Create a dask delayed object to gather the calibrations for the
    # instrument from each deployment
    cals = []
    for ind, row in deployments.iterrows():

        # Get the relevant request data
        deployNum, uid, deployStart, deployEnd = row["deploymentNumber"], row["uid"], row["deployStart"], row["deployEnd"]

        # Reformat the deployStart and deployEnd
        if pd.isna(deployEnd):
            deployEnd = datetime.datetime.now()

        # Construct the request parameters
        params = {
            "refdes": refdes,
            "uid": uid,
            "beginDT": deployStart.strftime("%Y-%m-%dT%H:%M:%S.000Z"),
            "endDT": deployEnd.strftime("%Y-%m-%dT%H:%M:%S.000Z")
        }

        # Request the data
        calInfo = get_api(URLS["cal"], params)

        # Reformat the data
        calDataFrame =  reformat_calInfo(calInfo, deployNum, uid)

        cals.append(calDataFrame)

    # Append individual calibrations into a single dataframe
    calibrations = pd.DataFrame()
    for cal in cals:
        calibrations = calibrations.append(cal, ignore_index=True)

    return calibrations


def get_calibrations_by_uid(uid, params=None):
    """Get calibrations for a specific instrument by its UID.

    Return the calibration coefficients for a given instrument
    in OOI for all of its deployments by its UID. The instrument
    may have been deployed at different locations (i.e. reference
    designators).

    Parameters
    ----------
    uid: (str)
        The unique id (UID) for the instrument for which to request
        vocab information.
    params: (dict), optional
        Optional parameters to pass to the calibration request to limit
        the data returned. Optional keywords are:
            beginDT: earliest datestring ("%Y-%m-%dT%H:%M:%S.000Z") to limit search
            endDT: latest datestring ("%Y-%m-%dT%H:%M:%S.000Z") to limit search

    Returns
    -------
    calibrations: (pandas.DataFrame)
        A table of the calibration information for the given instrument.
    """

    cal_url = URLS["cal"]
    if params is None:
        params = {
            "uid":uid
        }
    else:
        params["uid"] = uid

    # Make the calibration data request
    calInfo = get_api(cal_url, params=params)

    # Put the data into a pandas dataframe, sorted by calibration date and coefficient name
    columns = ["uid", "calCoef", "calDate", "value", "calFile"]
    calibrations = pd.DataFrame(columns=columns)
    for c in calInfo["calibration"]:
        for cc in c["calData"]:
            caldf = pd.DataFrame({
                "uid": cc["assetUid"],
                "calCoef": cc["eventName"],
                "calDate": convert_time(cc["eventStartTime"]),
                "value": cc["value"],
                "calFile": cc["dataSource"]
            })
            calibrations = pd.concat([calibrations, caldf], ignore_index=True)
    calibrations.sort_values(by=["calDate", "calCoef"], inplace=True)

    return calibrations

def get_deployments(refdes, deploy_num="-1", results=pd.DataFrame()):
    """Request deployment information for a reference designator.

    Get the deployment information for an instrument. Defaults to all
    deployments for a given instrument (reference designator) unless one is
    supplied.

    Parameters
    ----------
    refdes: (str)
        The reference designator for the instrument for which to request
        deployment information.
    deploy_num: (str)
        Optional to include a specific deployment number. Otherwise
        defaults to -1 which is all deployments.
    results: (pandas.DataFrame)
        Optional. Useful for recursive applications for gathering
        deployment information for multiple instruments.

    Returns
    -------
    results: (pandas.DataFrame)
        A table of the deployment information for the given instrument
        (reference designator) with deployment number, deployed water
        depth, latitude, longitude, start of deployment, end of deployment,
        and cruise IDs for the deployment and recovery.
    """
    # First, build the request
    array, node, instrument = refdes.split("-", 2)
    deploy_url = "/".join((URLS["deploy"], array, node, instrument,
                           deploy_num))

    # Next, get the deployments from the deploy url. The API returns a list
    # of dictionary objects with the deployment data.
    deployments = get_api(deploy_url)

    # Now, iterate over the deployment list and get the associated data for
    # each individual deployment
    while len(deployments) > 0:
        # Get a single deployment
        deployment = deployments.pop()

        # Process the dictionary data
        # Deployment Number
        deploymentNumber = deployment.get("deploymentNumber")

        # Location info
        location = deployment.get("location")
        depth = location["depth"]
        lat = location["latitude"]
        lon = location["longitude"]

        # Sensor info
        sensor = deployment.get("sensor")
        uid = sensor["uid"]
        assetId = sensor["assetId"]

        # Start and end times of the deployments
        startTime = convert_time(deployment.get("eventStartTime"))
        stopTime = convert_time(deployment.get("eventStopTime"))

        # Cruise IDs of the deployment and recover cruises
        deployCruiseInfo = deployment.get("deployCruiseInfo")
        recoverCruiseInfo = deployment.get("recoverCruiseInfo")
        if deployCruiseInfo is not None:
            deployID = deployCruiseInfo["uniqueCruiseIdentifier"]
        else:
            deployID = None
        if recoverCruiseInfo is not None:
            recoverID = recoverCruiseInfo["uniqueCruiseIdentifier"]
        else:
            recoverID = None

        # Put the data into a pandas dataframe
        data = np.array([[deploymentNumber, uid, assetId, lat, lon, depth,
                          startTime, stopTime, deployID, recoverID]])
        columns = ["deploymentNumber", "uid", "assetId",  "latitude",
                   "longitude", "depth", "deployStart", "deployEnd",
                   "deployCruise", "recoverCruise"]
        df = pd.DataFrame(data=data, columns=columns)

        # Generate the table results
        results = pd.concat([results, df])

    # Sort the deployments by deployment number and reset the index
    results = results.sort_values(by="deploymentNumber")
    results = results.reset_index(drop=True)

    return results


def get_datastreams(refdes):
        """Retrieve methods and data streams for a reference designator."""
        # Build the url
        array, node, instrument = refdes.split("-", 2)
        method_url = "/".join((URLS["data"], array, node, instrument))

        # Build a table linking the reference designators, methods, and data
        # streams
        stream_df = pd.DataFrame(columns=["refdes", "method", "stream"])
        methods = get_api(method_url)
        for method in methods:
            if "bad" in method:
                continue
            stream_url = "/".join((method_url, method))
            streams = get_api(stream_url)
            new_dict = {
                "refdes": refdes,
                "method": method,
                "stream": streams,
            }
            stream_df = pd.concat([stream_df, pd.DataFrame(new_dict)])

        # Expand so that each row of the dataframe is unique
        stream_df = stream_df.explode('stream').reset_index(drop=True)

        # Return the results
        return stream_df



def parse_metadata(metadata):
    """Parse metadata to dataframe.

    Parse the metadata dictionary for an instrument returned by OOI into
    a pandas dataframe.
    """
    # Put the two keys into separate dataframes
    metadata_times = pd.DataFrame(metadata["times"])
    metadata_parameters = pd.DataFrame(metadata["parameters"])

    # Merge the two into a single dataframe
    results = metadata_parameters.merge(metadata_times, left_on="stream",
                                        right_on="stream")
    results.drop_duplicates(inplace=True)

    # Return the results
    return results


def get_metadata(refdes):
    """Request metadata.

    Get the OOI Metadata for a specific instrument specified by its
    associated reference designator.

    Parameters
    ----------
    refdes: (str)
        OOINet standardized reference designator in the form of
        <array>-<node>-<instrument>.

    Returns
    -------
    results: (pandas.DataFrame)
        A dataframe with the relevant metadata of the given reference
        designator.
    """
    # First, construct the metadata request url
    array, node, instrument = refdes.split("-", 2)
    metadata_request_url = "/".join((URLS["data"], array, node,
                                     instrument, "metadata"))

    # Request the metadata
    metadata = get_api(metadata_request_url)

    # Parse the metadata
    metadata = parse_metadata(metadata)

    # Add in the reference designator
    metadata["refdes"] = refdes

    # Return the metadata
    return metadata


def get_parameter_data_levels(metadata):
    """Get parameters processing levels.

    Get the data levels (processing level) associated with the parameters
    for a given reference designator.

    Parameters
    ----------
    metadata: (pandas.DataFrame)
        The dataframe of metadata returned by get_metadata which contains
        the metadata for a given reference designator.

    Returns
    -------
    pid_dict: (dict)
        A dictionary with the data levels for each parameter id (Pid)
    """
    pdIds = np.unique(metadata["pdId"])
    pid_dict = {}
    for pid in pdIds:
        # Build the preload url
        preload_url = "/".join((URLS["preload"], pid.strip("PD")))
        # Query the preload data
        preload_data = get_api(preload_url)
        data_level = preload_data.get("data_level")
        # Update the results dictionary
        pid_dict.update({pid: data_level})

    return pid_dict

def get_annotations(refdes, **kwargs):
    """Retrieve data annotations for a given reference designator.

    Parameters
    ----------
    refdes: (str)
        The reference designator which to query OOINet for the
        associated annotation data.

    Kwargs
    ------
    method: (str)
        An OOINet method associated with the reference designator.
        Limits annotations for the given reference designator
        to only that method.
    stream: (str)
        An OOINet stream associated with the reference designator.
        Limits annotations for the given reference designator
        to only that stream.
    beginDT: (str)
        Limit the data request to only data after this date. Date
        should be formatted using OOINet epoch time (can be
        calculated using _)
    endDT: (str)
        Limit the data request to only data before this date.
    """
    # Need to build a parameters dictionary to pass to requests
    params = {"refdes": refdes}
    for key in kwargs:
        val = kwargs.get(key)
        # Convert datetimes to unix epoch
        if key == "beginDT":
            val = unix_epoch_time(val)
        elif key == "endDT":
            val = unix_epoch_time(val)
        else:
            pass
        params.update({key: val})

    # Get the annotations as a json and put into a dataframe
    anno_data = get_api(URLS["anno"], params=params)
    anno_data = pd.DataFrame(anno_data)

    # Convert the flags to QARTOD flags
    codes = {
        None: 0,
        'pass': 1,
        'not_evaluated': 2,
        'suspect': 3,
        'fail': 4,
        'not_operational': 9,
        'not_available': 9,
        'pending_ingest': 9
    }
    anno_data['qcFlag'] = anno_data['qcFlag'].map(codes).astype('category')

    return anno_data


def get_thredds_url(refdes, method, stream, goldCopy=False, **kwargs):
    """
    Return the url for the THREDDS server for the desired dataset(s).

    Parameters
    ----------
    refdes: (str)
        Reference designator for the instrument
    method: (str)
        The method (i.e. telemetered) for the given reference designator
    stream: (str)
        The stream associated with the reference designator and method
    goldCopy: (boolean: default False)
        Determine if should request from ooinet or should utilize the
        goldCopy THREDDS server.

    Kwargs
    ------
    beginDT: (str)
        Limit the data request to only data after this date.
    endDT: (str)
        Limit the data request to only data before this date.
    format: (str)
        e.g. "application/netcdf" (the default)
    include_provenance (str):
        'true' returns a text file with the provenance information
    include_annotations: (str)
        'true' returns a separate text file with annotations for the data
        within the given date range

    Returns
    -------
    thredds_url: (str)
        A url to the OOI Thredds server which contains the desired datasets
    """
    # Build the data request url
    if goldCopy is True:
        goldCopy_url = URLS["goldCopy"] + "-".join((refdes, method, stream)) + "/catalog.html"
        return goldCopy_url
    else:
        array, node, instrument = refdes.split("-", 2)
        data_request_url = "/".join((URLS["data"], array, node,
                                 instrument, method, stream))

    # Ensure proper datetime format for the request
    if 'beginDT' in kwargs.keys():
        kwargs['beginDT'] = pd.to_datetime(kwargs['beginDT']).strftime(
            '%Y-%m-%dT%H:%M:%S.%fZ')
    if 'endDT' in kwargs.keys():
        kwargs['endDT'] = pd.to_datetime(kwargs['endDT']).strftime(
            '%Y-%m-%dT%H:%M:%S.%fZ')

    # Build the query
    if len(kwargs) == 0:
        kwargs = {"require_deployment":True}
    else:
        kwargs.update({
            "require_deployment":True
        })
    params = kwargs

    # Request the data
    print("Waiting for request to process")
    # Get the urls
    urls = get_api(data_request_url, params=params)
    # Check the status of the dataset preparation
    status_url = [url for url in urls["allURLs"] if re.match(r'.*async_results.*', url)][0]
    status_url = status_url + "/status.txt"
    status = SESSION.get(status_url)
    dt = 0
    while status.status_code != requests.codes.ok:
        time.sleep(2)
        dt += 2
        if dt%2 == 0 and dt%5 == 0:
            print("Waiting for request to process")
        status = SESSION.get(status_url)

    # The asynchronous data request is contained in the 'allURLs' key,
    # in which we want to find the url to the thredds server
    for d in urls['allURLs']:
        if 'thredds' in d:
            thredds_url = d

    return thredds_url


def get_thredds_catalog(thredds_url):
    """
    Get the dataset catalog for the requested data stream.

    Parameters
    ----------
    thredds_url (str): the THREDDS server url for the
        requested data stream

    Returns
    -------
    catalog (list): the THREDDS catalog of datasets for
        the requested data stream
    """
    # Parse out the dataset_id from the thredds url
    page = SESSION.get(thredds_url).text
    soup = BeautifulSoup(page, "html.parser")
    pattern = re.compile('.*\\.nc$')
    catalog = sorted([node.get('href') for node in soup.find_all('a', text=pattern)])

    # Return the catalog
    return catalog


def download_netCDF_files(catalog, goldCopy=False, saveDir=None, verbose=True):
    """Download netCDF files for given netCDF datasets.

    Downloads the netCDF files returned by parse_catalog. If no path is
    specified for the save directory, will download the files to the
    current working directory.

    Parameters
    ----------
    catalog: (list)
        The netCDF catalog which to download
    goldCopy: (boolean)
        If you should use the gold copy thredds server (which is static)
    saveDir: (str)
        The path to the directory where to download the netCDF files. If no path
        is specified will download to your current directory
    verbose: (boolean)
        If True, show status of download.
    """
    # Step 1 - parse the catalog to get the correct web address to download the files
    if goldCopy is True:
        fileServer = URLS["goldCopy_fileServer"]
    else:
        fileServer = URLS["fileServer"]
    netCDF_files = [re.sub("catalog.html\?dataset=", fileServer, file) for file in catalog]

    # Step 2 - Specify and make the relevant save directory
    setup_download_dir(saveDir)

    # Step 3 - Initialize workers
    # 3A. Identify number availabe cores and divide in half
    cores = os.cpu_count()
    workers = int(cores / 2)
    # 3B. Setup a Queue
    queue = Queue()
    # 3C. Start workers
    for worker in range(workers):
        worker = DownloadWorker(queue)
        worker.daemon = True
        worker.start()

    # Step 4. Execute the processes and download the files
    print("----- Downloading files -----")
    for file in netCDF_files:
        queue.put((saveDir, file))
    # Execute
    queue.join()


def clean_catalog(catalog, stream, deployments=None):
    """Clean up the THREDDS catalog of unwanted datasets"""
    # Next, check that the netCDF datasets are not empty by getting the timestamps in the
    # datasets and checking if they are 
    datasets = []
    for dset in catalog:
        # Get the timestamps
        timestamps = dset.split("_")[-1].replace(".nc","").split("-")
        t1, t2 = timestamps
        # Check if the timestamps are equal
        if t1 == t2:
            pass
        else:
            datasets.append(dset)
            
    # Next, check if you want to filter for certain deployments
    if deployments is not None:
        # Check that the deployments are a list
        deployments = list(deployments["deploymentNumber"].astype(int))
        datasets = []
        for dset in catalog:
            dep = re.findall("deployment[\d]{4}", dset)[0]
            depNum = int(dep[-4:])
            if depNum not in deployments:
                pass
            else:
                datasets.append(dset)        
        
    # Finally, determine if the dataset is either for the given instrument
    # or an ancillary instrument which supplies and input variable
    catalog = datasets
    datasets = []
    ancillary = []
    for dset in catalog:
        check = dset.split("/")[-1]
        if stream in check:
            datasets.append(dset)
        else:
            ancillary.append(dset)
                       
    return datasets, ancillary
