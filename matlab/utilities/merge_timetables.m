function [merged, fail] = merge_timetables(timetables)
% MERGE_TIMETABLES Merge multiple timetables approximating the xarray
% compat='override' methodology used in the corresponding Python code base.
%
% cwingard 2026-01-06 -- initially developed using Claude AI 4.5 with
%   revision and refinement based on testing

if isempty(timetables)
    error('merge_timetables:EmptyInput', 'Input cell array is empty');
end %if

% remove any elements that are not timetables
m = ~cellfun(@istimetable, timetables);
timetables(m) = [];
if isempty(timetables)
    error('merge_timetables:ValueError', 'None of the objects in the array were timetables');
end %if
clear m

% Start with the first timetable
tic
merged = timetables{1};
ntables = numel(timetables);
fail = 0;

% set up reporting metrics
percentage = 1;
backspaces = '';
percent_step = 100 / (ntables - 1);

% Merge each subsequent timetable
for i = 2:ntables
    try
        merged = merge_two_timetables(merged, timetables{i});
    catch
        fail = fail + 1;
    end %try
    % Print percentage progress
    percentage = percentage + percent_step;
    perc_str = sprintf('... percent merged: %3.1f', percentage);
    fprintf([backspaces, perc_str]);
    backspaces = repmat(sprintf('\b'), 1, length(perc_str));
end %for
fprintf("\n")
clear i percentage backspaces percent_step perc_str

% check the fail count
if fail == ntables - 1
    % if fail equals the number of table objects, then most likely the first table is bad
    fprintf('First merging attempt failed, retrying ...\n')
    timetables{1} = [];  % drop the first table
    merged = timetables{1};  % get set to try again

    ntables = numel(timetables);
    percentage = 1;
    backspaces = '';
    percent_step = 100 / (ntables - 1);
    fail = 1;
    for i = 2:ntables
        try
            merged = merge_two_timetables(merged, timetables{i});
        catch
            fail = fail + 1;
        end %try
        % Print percentage progress
        percentage = percentage + percent_step;
        perc_str = sprintf('... percent merged: %3.1f', percentage);
        fprintf([backspaces, perc_str]);
        backspaces = repmat(sprintf('\b'), 1, length(perc_str));
    end %for
    fprintf("\n")
    clear i percentage backspaces percent_step perc_str
    
    if fail == ntables - 1
        error('merge_timetables:ValueError', 'Unable to merge the data, please review source data for possible errors.');
    end %if
end %if
clear ntables
toc
end %function


function merged = merge_two_timetables(tt1, tt2)
% MERGE_TWO_TIMETABLES Merges two timetables together, adding missing
% variables using datatype specific fill values
%
% cwingard 2026-01-06 -- initially developed using Claude AI 4.5 with
%   revision and refinement based on testing

% Get variable names
vars1 = tt1.Properties.VariableNames;
vars2 = tt2.Properties.VariableNames;

% Find overlapping variables
commonVars = intersect(vars1, vars2);

% Standardize dimensions for common variables (take larger size, pad with fill value)
for v = 1:length(commonVars)
    varName = commonVars{v};
    data1 = tt1.(varName);
    data2 = tt2.(varName);
    
    [~, cols1] = size(data1);
    [~, cols2] = size(data2);
    
    if cols1 ~= cols2
        % Use the larger dimension
        targetCols = max(cols1, cols2);
        
        % Pad tt1 if needed
        if cols1 < targetCols
            fillVal = get_fill_value(data1);
            padding = repmat(fillVal, size(data1, 1), targetCols - cols1);
            tt1.(varName) = [data1, padding];
        end
        
        % Pad tt2 if needed
        if cols2 < targetCols
            fillVal = get_fill_value(data2);
            padding = repmat(fillVal, size(data2, 1), targetCols - cols2);
            tt2.(varName) = [data2, padding];
        end
    end
end
    
% initialize the merged table based on time and deployment number
merged = cat(1, tt1(:, 'deployment'), tt2(:, 'deployment'));

% Add variables from tt1
[~, idx] = ismember(tt1(:, "deployment"), merged(:, "deployment"));
for v = 1:length(vars1)
    varName = vars1{v};
    if varName == "deployment"
        continue
    end %if
    varData = tt1.(varName);
    
    % Initialize with appropriate fill value based on data type
    [~, numCols] = size(varData);
    mergedData = initialize_with_fill_value(varData, height(merged), numCols);
    
    % Fill in values from tt1
    mergedData(idx, :) = varData;
    merged.(varName) = mergedData;
end

% update the metadata attributes in merged
merged = copy_timetable_metadata(tt1, merged);

% Add variables from tt2
[~, idx] = ismember(tt2(:, "deployment"), merged(:, "deployment"));
for v = 1:length(vars2)
    varName = vars2{v};
    if varName == "deployment"
        continue
    end %if
    varData = tt2.(varName);
    [~, numCols] = size(varData);

    if ~ismember(varName, vars1)
        % New variable - initialize with appropriate fill value
        mergedData = initialize_with_fill_value(varData, height(merged), numCols);
    else
        % Existing variable
        mergedData = merged.(varName);
    end
    
    % add the values from tt2
    mergedData(idx, :) = varData;    
    merged.(varName) = mergedData;
end

% update the metadata attributes in merged
merged = copy_timetable_metadata(tt2, merged);
end %function


function mergedData = initialize_with_fill_value(varData, numRows, numCols)
% INITIALIZE_WITH_FILL_VALUE Initialize array with appropriate fill value based
% on data type
%
% cwingard 2026-01-06 -- initially developed using Claude AI 4.5 with
%   revision and refinement based on testing

if isnumeric(varData)
    % Use -Inf for all numeric types (float and integer)
    mergedData = repmat(-Inf, numRows, numCols);
elseif ischar(varData)
    % Empty string for char arrays
    mergedData = repmat('', numRows, numCols);
elseif isstring(varData)
    % Empty string for string arrays
    mergedData = strings(numRows, numCols);
elseif isdatetime(varData)
    % NaT for datetime
    mergedData = repmat(NaT, numRows, numCols);
else
    % Generic missing for other types
    mergedData = repmat(varData(1,:), numRows, 1);
    mergedData(:) = missing;
end %if
end %function


function fillVal = get_fill_value(varData)
% GET_FILL_VALUE Get a single fill value appropriate for the data type
%
% cwingard 2026-01-06
dataType = class(varData);    
if isnumeric(varData)
    % Use -Inf for all numeric types (float and integer)
    fillVal = cast(-Inf, dataType);
elseif ischar(varData)
    % Empty string for char
    fillVal = '';
elseif isstring(varData)
    % Empty string for string
    fillVal = "";
elseif isdatetime(varData)
    % NaT for datetime
    fillVal = NaT;
else
    % Generic missing for other types
    fillVal = missing;
end
end


function dest_timetable = copy_timetable_metadata(source_timetable, dest_timetable)
% COPY_TIMETABLE_METADATA Copy metadata from source timetable to destination timetable
%
% cwingard 2026-01-06 -- initially developed using Claude AI 4.5 with
%   revision and refinement based on testing

% Copy timetable-level properties
if ~isempty(source_timetable.Properties.Description)
    dest_timetable.Properties.Description = source_timetable.Properties.Description;
end

if ~isempty(source_timetable.Properties.UserData)
    dest_timetable.Properties.UserData = source_timetable.Properties.UserData;
end

% Copy DimensionNames if they exist
if ~isempty(source_timetable.Properties.DimensionNames)
    dest_timetable.Properties.DimensionNames = source_timetable.Properties.DimensionNames;
end

% Copy variable-specific metadata for matching variable names
source_vars = source_timetable.Properties.VariableNames;
dest_vars = dest_timetable.Properties.VariableNames;

% Find common variables
[common_vars, source_idx, dest_idx] = intersect(source_vars, dest_vars, 'stable');

% Copy VariableUnits for common variables
if ~isempty(source_timetable.Properties.VariableUnits)
    for i = 1:length(common_vars)
        dest_timetable.Properties.VariableUnits{dest_idx(i)} = ...
            source_timetable.Properties.VariableUnits{source_idx(i)};
    end
end

% Copy VariableDescriptions for common variables
if ~isempty(source_timetable.Properties.VariableDescriptions)
    for i = 1:length(common_vars)
        dest_timetable.Properties.VariableDescriptions{dest_idx(i)} = ...
            source_timetable.Properties.VariableDescriptions{source_idx(i)};
    end
end

% Copy VariableContinuity for common variables
if ~isempty(source_timetable.Properties.VariableContinuity)
    for i = 1:length(common_vars)
        dest_timetable.Properties.VariableContinuity{dest_idx(i)} = ...
            source_timetable.Properties.VariableContinuity{source_idx(i)};
    end
end

% Copy custom properties (if any exist beyond standard ones)
source_props = properties(source_timetable.Properties);
dest_props = properties(dest_timetable.Properties);

standard_props = {'Description', 'UserData', 'DimensionNames', ...
                  'VariableNames', 'VariableTypes', 'VariableUnits', ...
                  'VariableDescriptions', 'VariableContinuity', ...
                  'RowTimes', 'StartTime', 'SampleRate', 'TimeStep', ...
                  'CustomProperties', 'Events'};

custom_props = setdiff(source_props, standard_props);

for i = 1:length(custom_props)
    prop_name = custom_props{i};
    if ismember(prop_name, dest_props)
        try
            dest_timetable.Properties.(prop_name) = source_timetable.Properties.(prop_name);
        catch ME
            warning('Could not copy property %s: %s', prop_name, ME.message);
        end
    end
end
end