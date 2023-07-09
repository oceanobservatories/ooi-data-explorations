function load_kdata()

def load_kdata(site, node, sensor, method, stream, tag='*.nc', use_dask=False):
    """
    Download data from the JupyterHub kdata directories, using the reference
    designator parameters to select the catalog of interest and the regex tag
    to select the NetCDF files of interest. In most cases, the default tag
    can be used, however for instruments that require data from a co-located
    sensor, a more detailed regex tag will be required to ensure that only data
    files from the instrument of interest are loaded.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :param method: Delivery method for the data (either telemetered,
        recovered_host or recovered_inst)
    :param stream: Stream name that contains the data of interest
    :param tag: regex pattern to select the NetCDF files to download
    :param use_dask: Boolean flag indicating whether to load the data using
        dask arrays (default=False)
    :return data: All the data, combined into a single dataset
    """
    # download the data from the Gold Copy THREDDS server
    dataset_id = '-'.join([site, node, sensor, method, stream])
    data = kdata_collect(dataset_id, tag, use_dask)
    return data


def kdata_collect(dataset_id, tag='*.nc', use_dask=False):
    """
    Use a regex tag combined with the dataset ID to collect data from the OOI
    JupyterHub kdata directory. The collected data is gathered into a xarray
    dataset for further processing.

    :param dataset_id: dataset ID as a string
    :param tag: regex tag to use in discriminating the data files, so we only
        collect the data files of interest
    :param use_dask: Boolean flag indicating whether to load the data using
        dask arrays (default=False)
    :return gc: the collected Gold Copy data as a xarray dataset
    """
    # construct the kdata directory path with the dataset ID
    kdata = os.path.abspath(os.path.join(os.path.expanduser('~'), 'ooi/kdata'))
    kdata = os.path.abspath(os.path.join(kdata, dataset_id))

    # Create a list of the files from the request above using a simple regex as a tag to discriminate the files
    files = glob.glob(kdata + '/' + tag)

    # Process the data files found above and concatenate them into a single list
    print('Downloading %d data file(s) from the local kdata directory' % len(files))
    if len(files) < 4:
        # just 1 to 3 files, download sequentially
        frames = [process_file(file, gc='KDATA', use_dask=use_dask) for file in tqdm(files, desc='Loading and '
                                                                                                 'Processing Data '
                                                                                                 'Files')]
    else:
        # multiple files, use multithreading to download concurrently
        part_files = partial(process_file, gc='KDATA', use_dask=use_dask)
        with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
            frames = list(tqdm(executor.map(part_files, files), total=len(files),
                               desc='Loading and Processing Data Files', file=sys.stdout))

    if not frames:
        message = "No data files were loaded from the JupyterHub kdata directory."
        warnings.warn(message)
        return None

    # merge the data frames into a single data set
    data = merge_frames(frames)

    return data
