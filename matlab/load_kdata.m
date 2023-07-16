function data = load_kdata(site, node, sensor, method, stream, tag)
% LOAD_KDATA Download data from the OOI kdata directory
%
% Download data from the local kdata directory available on the OOI JupyterHub
% (the files in this directory are served by the OOI Gold Copy THREDDS catalog)
% using reference designator parameters to select the folder of interest and a
% glob pattern to select the NetCDF files of interest.
% 
% INPUTS:
%
%   site -- Site designator, extracted from the first part of the
%       reference designator
%   node -- Node designator, extracted from the second part of the
%       reference designator
%   sensor -- Sensor designator, extracted from the third and fourth part
%       of the reference designator
%   method -- Delivery method for the data (either telemetered,
%       recovered_host or recovered_inst)
%   stream -- Stream name that contains the data of interest
%   tag --  glob pattern to select the NetCDF files to download. This should
%       NOT be too broad (e.g., "*.nc") or you will be downloading all of the
%       data files in a catalog! Instead, refine the glob pattern (e.g., add
%       deployment number and the instrument class name, 
%       "deployment0003*OPTAA*.nc") to limit the data downloaded.
%
% OUTPUTS:
%
%   data -- All the data for the particular dataset of interest, combined into
%       a single timetable
%
% C. Wingard, 2023-07-10

% load the common file processing utilities
ph = process_files;

% set up the name of kdata directory (full path) and the dataset ID
base_dir = "/home/jovyan/ooi/kdata/";
if ~isfolder(base_dir)
    error("This function can only be used on the OOI JupyterHub")
end %if
dataset_id = join([upper(string(site)), upper(string(node)), upper(string(sensor)), lower(string(method)), lower(string(stream))], "-");

% determine the number of available files to download based on the dataset ID
% and the regex tag
data_dir = join([base_dir dataset_id], "");
files = dir(join([data_dir "/" tag], ""));
nfiles = numel(files);      % determine number of files to download
data = cell(nfiles, 1);     % pre-allocate a cell array for downloaded results

% check to see if the user has the parallel processing toolbox and there are
% enough files to justify its use (buildup and teardown of the parallel
% processing utilities costs time, so we only want to use it when it benefits
% us)
fprintf('Loading %d files from the local kdata directory ...\n', nfiles)
tic
% use a traditional for loop to download the files, adding them to the 
% pre-allocated cell array
for i = 1:nfiles
    nc_file = join([data_dir "/" files(i).name], "");
    data{i} = ph.process_file(nc_file);
end %for
toc
clear url files nfiles dap_url repetition file ext nc_url

% finally, concatenate the timetables together returning the results
data = ph.merge_frames(data);
end %function
