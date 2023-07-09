function load_gc_thredds(site, node, sensor, method, stream, tag)
% Download data from the OOI Gold Copy THREDDS catalog, using the reference
% designator parameters to select the catalog of interest and the regex tag
% to select the NetCDF files of interest. In most cases, the default tag
% can be used, however for instruments that require data from a co-located
% sensor, a more detailed regex tag will be required to ensure that only data
% files from the instrument of interest are loaded.
% 
% :param site: Site designator, extracted from the first part of the
%     reference designator
% :param node: Node designator, extracted from the second part of the
%     reference designator
% :param sensor: Sensor designator, extracted from the third and fourth part
%     of the reference designator
% :param method: Delivery method for the data (either telemetered,
%     recovered_host or recovered_inst)
% :param stream: Stream name that contains the data of interest
% :param tag: regex pattern to select the NetCDF files to download
% :return data: All the data, combined into a single dataset

% load the common file processing utilities
ph = process_files;

% set up the GC THREDDS catalog URL, the dataset ID and a default regex tag it
% it wasn't set
base_url = "https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/";
dataset_id = join([upper(site), upper(node), upper(sensor), lower(method), lower(stream)], "-");
if ~input('tag', 'var')
    tag = ".*\.nc$";
end %if

% determine the number of available files to download based on the dataset ID
% and the regex tag
url = join([base_url dataset_id "/catalog.html"], "");
files = ph.list_files(url, tag);

% change the URL to the file server and start downloading data files
file_url = 'https://thredds.dataexplorer.oceanobservatories.org/thredds/fileServer/';
data = cell(size(files, 1), 1);
p = parpool("Threads");
for repetition = 1:numel(files)
    tmpName = [tempname, '.nc'];
    nc_url = join([file_url, files(repetition)], '');
    websave(tmpName, nc_url);
    data{repetition} = ph.process_file(tmpName);
    delete(tmpName)
end %if
delete(p)
clear file_url p repetition tmpName

end %functions

% def gc_collect(dataset_id, tag='.*\\.nc$', use_dask=False):
%     """
%     Use a regex tag combined with the dataset ID to collect data from the OOI
%     Gold Copy THREDDS catalog. The collected data is gathered into a xarray
%     dataset for further processing.
% 
%     :param dataset_id: dataset ID as a string
%     :param tag: regex tag to use in discriminating the data files, so we only
%         collect the data files of interest
%     :param use_dask: Boolean flag indicating whether to load the data using
%         dask arrays (default=False)
%     :return gc: the collected Gold Copy data as a xarray dataset
%     """
%     # construct the THREDDS catalog URL based on the dataset ID
%     gc_url = 'http://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/'
%     url = gc_url + dataset_id
% 
%     # Create a list of the files from the request above using a simple regex as a tag to discriminate the files
%     files = list_files(url, tag)
% 
%     # Process the data files found above and concatenate them into a single list
%     print('Downloading %d data file(s) from the OOI Gold Copy THREDSS catalog' % len(files))
%     if len(files) < 4:
%         # just 1 to 3 files, download sequentially
%         frames = [process_file(file, gc='GC', use_dask=use_dask) for file in tqdm(files, desc='Downloading and '
%                                                                                               'Processing Data '
%                                                                                               'Files')]
%     else:
%         # multiple files, use multithreading to download concurrently
%         part_files = partial(process_file, gc='GC', use_dask=use_dask)
%         with ThreadPoolExecutor(max_workers=N_CORES) as executor:
%             frames = list(tqdm(executor.map(part_files, files), total=len(files),
%                                desc='Downloading and Processing Data Files', file=sys.stdout))
% 
%     if not frames:
%         message = "No data files were downloaded from the Gold Copy THREDDS server."
%         warnings.warn(message)
%         return None
% 
%     # merge the data frames into a single data set
%     data = merge_frames(frames)
% 
%     return data
