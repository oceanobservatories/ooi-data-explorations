function fh = process_files()
% PROCESS_FILES Function handles to common files used to load NetCDF data
%
% cwingard 2023-07-09
fh.list_files = @list_files;
fh.process_file = @process_file;
fh.merge_frames = @merge_frames;
end %function

function nc_files = list_files(url, tag)
% LIST_FILES create a list of the NetCDF data files in a THREDDS catalog
% Inputs are the URL for the THREDDS catalog of interest and a regex tag
% used to subselect only the NetCDF files of interest
%
% cwingard 2023-07-09

% use webread to scrape the catalog
catalog = webread(url);

% pull a list of netCDF files out of the catalog
nc_all = regexp(catalog, '<a href=''([^>]+.nc)''>', 'tokens');

% cell elements are themselves cells, convert to cell array of character vectors
nc_all = [nc_all{:}]';

% drop 'catalog.html?dataset=' from URL
nc_all = cellfun(@(x) x(22:end), nc_all, 'UniformOutput', 0); %

% use the regex tag to select only files of interest
nc_all = cellfun(@(x) regexp(x, tag, 'match'), nc_all, 'UniformOutput', 0);
nc_all(cellfun('isempty', nc_all)) = [];  % remove empty matches

% convert the list to a string array
nc_files = string(nc_all(:));

% check to see if we found any files to download, if not throw an error
if size(nc_files, 1) == 0
    error("Unable to find any files to download using %s at %s", [tag, url])
end %if
end %function

function data = process_file(filename)
% PROCESS_FILE Load local NetCDF into a timetable
%
% cwingard 2023-07-09

% Get information about the NetCDF data file
file_info = h5info(filename);

% load the data into a timetable
data = nc_reader(filename);

% remove some of the variables that are not used or are better served elsewhere
m = ismember(data.Properties.VariableNames, {'id', 'provenance', ...
    'dcl_controller_timestamp', 'driver_timestamp', 'ingestion_timestamp', ...
    'port_timestamp', 'preferred_timestamp', 'station', 'z'});
data = removevars(data, m);

% sort the data, first by the deployment number (if more than one, and then by
% the time
data = sortrows(data, {'deployment', 'Time'});

% add global attributes to the timetable (if this data came from OOI)
m = reshape(strcmp({file_info.Attributes.Name}, 'publisher_name'), size(file_info.Attributes));
if sum(m) == 1
    if strcmp(file_info.Attributes(m).Value, 'Ocean Observatories Initiative')
        % assign the NetCDF title attribute to the timetable description
        m = reshape(strcmp({file_info.Attributes.Name}, 'title'), size(file_info.Attributes));
        data.Properties.Description = file_info.Attributes(m).Value;
        
        % create custom properties for the timetable, pulling different global 
        % attributes out of the NetCDF file and assigning them
        prop_names = {'subsite', 'node', 'sensor', 'collection_method', ...
            'stream', 'date_created', 'date_modified', 'requestUUID', ...
            'Manufacturer', 'ModelNumber', 'Description', 'SerialNumber', ...
            'AssetUniqueID', 'AssetManagementRecordLastModified', ...
            'publisher_name', 'acknowledgement'};
        prop_types = {'table', 'table', 'table', 'table', 'table', 'table', ...
            'table', 'table', 'table', 'table', 'table', 'table', 'table', ...
            'table', 'table', 'table'};
        data = addprop(data, prop_names, prop_types);    
        for i = 1:size(prop_types, 2)
            m = reshape(strcmp({file_info.Attributes.Name}, prop_names{i}), ...
                size(file_info.Attributes));
            data.Properties.CustomProperties.(prop_names{i}) = file_info.Attributes(m).Value;
        end %for
        clear prop_names prop_types i
    end %if
end %if
clear file_info m

end %function

% def merge_frames(frames):
%     """
%     Merge the multiple data frames downloaded from the M2M system or the Gold
%     Copy THREDDS server into a single xarray data set. Keep track of how many
%     frames fail to merge.
% 
%     :param frames: The data frames to concatenate/merge into a single data set
%     :return data: The final, merged data set
%     """
%     # merge the list of processed data frames into a single data set
%     nframes = len(frames)
%     bad_frames = 0
%     if nframes > 1:
%         try:
%             # first try to just concatenate all the frames; this usually works, but not always
%             data = xr.concat(frames, dim='time')
%         except ValueError:
%             # try merging the frames one-by-one into a single data set
%             data, fail = _frame_merger(frames[0], frames)
% 
%             # if all files failed that would suggest the first file is the problem.
%             # try the merge again, resetting the starting frame to skip the first one.
%             if nframes - fail == 1:
%                 try:
%                     bad_frames += 1
%                     data = xr.concat(frames[1:], dim='time')
%                 except ValueError:
%                     # this data set has issues! try merging one more time, frame by frame
%                     data, fail = _frame_merger(frames[1], frames[1:])
%                     bad_frames += fail
% 
%                     # if we still can't merge the frames, then there probably is something more fundamentally wrong,
%                     # and trying to account for it here is not going to be possible
%                     if nframes - 1 - fail == 1:
%                         message = f"Unable to merge the {nframes} files downloaded from the Gold Copy THREDDS server."
%                         warnings.warn(message)
%                         return None
%             else:
%                 bad_frames += fail
%     else:
%         # there is just the one
%         data = frames[0]
% 
%     if bad_frames > 0:
%         message = "{} of the {} downloaded files failed to merge.".format(bad_frames, nframes)
%         warnings.warn(message)
% 
%     data = data.sortby(['deployment', 'time'])
%     _, index = np.unique(data['time'], return_index=True)
%     data = data.isel(time=index)
%     data.attrs['time_coverage_start'] = ('%sZ' % data.time.min().values)
%     data.attrs['time_coverage_end'] = ('%sZ' % data.time.max().values)
%     data.attrs['time_coverage_resolution'] = ('P%.2fS' % (np.mean(data.time.diff('time').values).astype(float) / 1e9))
% 
%     return data
% 
% 
% def _frame_merger(data, frames):
%     """
%     Internal method used by merge_frames to enumerate through the frames,
%     trying to concatenate/merge the data frames together into a single
%     data set.
% 
%     :param data: initial data frame to concatenate/merge with the other frames
%     :param frames: additional frames to add on to the initial data frame
%     :return data: the final concatenated/merged data set
%     :return fail: a count of the number of files that failed
%     """
%     fail = 0
%     for idx, frame in enumerate(frames[1:], start=2):
%         try:
%             # concatenation handles 99% of the cases
%             with dask.config.set(**{'array.slicing.split_large_chunks': False}):
%                 data = xr.concat([data, frame], dim='time')
%         except (ValueError, NotImplementedError):
%             try:
%                 # try merging the data, usually one of the data files is missing a variable from a co-located
%                 # sensor that the system was unable to find
%                 _, index = np.unique(data['time'], return_index=True)
%                 data = data.isel(time=index)
%                 with dask.config.set(**{'array.slicing.split_large_chunks': False}):
%                     data = data.merge(frame, compat='override')
%             except (ValueError, NotImplementedError):
%                 # something is just not right with this data file
%                 fail += 1
% 
%     return data, fail


