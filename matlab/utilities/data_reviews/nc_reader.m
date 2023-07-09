function t = nc_reader(filename)
% NC_READER Read a NetCDF file, extracting all variables (with a few
% exceptions), saving the results in a timetable.
% 
% Code based on example available from Mathworks at:
% https://www.mathworks.com/help/parallel-computing/examples/process-big-data-in-the-cloud.html

% Get information about the NetCDF data file
fileInfo = h5info(filename);

% Extract the global and variable level attributes -- note, Matlab doesn't
% really support these very well, so their utility is limited.
%gAttributes = struct2table(fileInfo.Attributes);
vAttributes = {fileInfo.Datasets.Attributes};

% Extract the variable names
varNames = string({fileInfo.Datasets.Name});

% test for presence of a variable called time
i = 1; test = 0;
while test == 0 && i <= numel(varNames)
    test = strcmp('time', varNames{i});
    i = i + 1;
end %while
if ~test
    error('The NetCDF file specified does not include the variable ''time''')
end %if
clear i test

% Create the datetime axis from the time variable (ERDDAP uses
% 1970, while the OOI-created NetCDF files use 1900 as their pivot years).
nc_time = h5read(filename, '/time');   % obtain the time record
test = nc_time(1) / 60 / 60 / 24 + datenum(1970, 1, 1, 0, 0, 0);
if test > now
    dt = datetime(1900, 1, 1, 'Format', 'dd-MMM-yyyy HH:mm:ss', ...
        'TimeZone', 'UTC') + seconds(nc_time);
else 
    dt = datetime(1970, 1, 1, 'Format', 'dd-MMM-yyyy HH:mm:ss', ...
        'TimeZone', 'UTC') + seconds(nc_time);
end %if
rowlength = length(dt);
clear test nc_time

% Create an empty timetable using the datetime axis
t = timetable('RowTimes', dt);

% Populate the timetable with the variable data
for k = 1:numel(varNames)
    % skip adding time (already added as the RowTime) and obs
    if strcmp(varNames{k}, 'time') || strcmp(varNames{k}, 'obs')
        continue
    end %if
    % read the variable from the NetCDF file
    data = squeeze(h5read(filename, "/" + varNames{k}));
    % pull out the variable units and comment attributes
    if ~isempty(vAttributes{k})
        attr = struct2table(vAttributes{k});
        units = {''}; descr = {''};
        for j = 1:height(attr)
            if strcmp(attr.Name(j), 'units')
                units = attr.Value(j);
            end %if
            if strcmp(attr.Name(j), 'comment')
                descr = attr.Value(j);
            end %if
        end %for
    end %if
    % add the variable and attributes to the time table
    [r, c] = size(data);    % check the dimensions 
    if r == rowlength
        % if the number of rows == the number of RowTimes, add the variable
        % without modification.
        if strcmp(fileInfo.Datasets(k).Datatype.Class, 'H5T_STRING')
            t = addvars(t, cell2mat(data), 'NewVariableNames', varNames{k});
        else
            t = addvars(t, data, 'NewVariableNames', varNames{k});
        end %if
        t.Properties.VariableUnits{varNames{k}} = units{:};
        t.Properties.VariableDescriptions{varNames{k}} = descr{:};
    elseif c == rowlength
        % if the number of columns equals the RowTimes, rotate the variable
        % before adding it so the row length matches the RowTimes
        if strcmp(fileInfo.Datasets(k).Datatype.Class, 'H5T_STRING')
            t = addvars(t, cell2mat(data'), 'NewVariableNames', varNames{k});
        else
            t = addvars(t, data', 'NewVariableNames', varNames{k});
        end %if
        t.Properties.VariableUnits{varNames{k}} = units{:};
        t.Properties.VariableDescriptions{varNames{k}} = descr{:};
    elseif r == 1 && c == 1
        % this is a scalar variable, and it needs to replicated out to the
        % RowTimes dimension before it can be added.
        t = addvars(t, repmat(data, rowlength, 1), 'NewVariableNames', varNames{k});
        t.Properties.VariableUnits{varNames{k}} = units{:};
        t.Properties.VariableDescriptions{varNames{k}} = descr{:};        
    else
        % this is something weird, ignore it for now.
    end %if
end %for
clear dt rowlength k data r c
end %function
