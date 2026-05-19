function ds = add_annotation_qc_flags(ds, annotations)
% ADD_ANNOTATION_QC_FLAGS Adds the annotation QC flags to a timetable as
% new variable(s).
%
% cwingard 2026-01-06 -- adapted from the Python function of the same name

% Convert the flags to QARTOD flags
codes = containers.Map({'note', 'pass', 'not_evaluated', 'suspect', 'fail', ...
                       'not_operational', 'not_available', 'pending_ingest'}, ...
                       [0, 1, 2, 3, 4, 9, 9, 9]);

% Map qcFlag column using the codes dictionary
for i = 1:height(annotations)
    if ismissing(annotations.qcFlag(i)) || isempty(annotations.qcFlag{i})
        annotations.qcFlag(i) = {0};  % No flag set
    elseif isKey(codes, annotations.qcFlag{i})
        annotations.qcFlag(i) = {codes(annotations.qcFlag{i})};
    end
end
annotations.qcFlag = str2double(annotations.qcFlag);

% Filter only for annotations which apply to the dataset
stream = ds.Properties.CustomProperties.stream;
streamMask = cellfun(@(x) isempty(x) || strcmp(x, stream), annotations.stream);
annotations = annotations(streamMask, :);

% Explode the annotations so each parameter is hit for each annotation
% In MATLAB, we need to manually expand rows with cell arrays
expandedRows = {};
for i = 1:height(annotations)
    newRow = annotations(i, :);
    params = annotations.parameters(i);
    if ismissing(params)
        newRow.parameters = "";
        expandedRows{end+1} = newRow; %#ok<AGROW>
        continue
    end %if
    params = num2cell(jsondecode(params));
    for j = 1:length(params)
        newRow.parameters = params(j);
        expandedRows{end+1} = newRow; %#ok<AGROW>
    end
end
annotations = vertcat(expandedRows{:});

% Get the unique parameters and their associated variable name
streamAnnos = containers.Map();
uniqueParams = unique([annotations.parameters]);
for i = 1:length(uniqueParams)
    pid = str2double(uniqueParams(i));
    if isnan(pid)
        paramName = 'rollup';
    else
        paramInfo = get_parameter_information(char(uniqueParams(i)));
        paramName = paramInfo.netcdf_name;
    end
    streamAnnos(paramName) = pid;
end

% Next, get the flags associated with each parameter or all parameters
flagsDict = containers.Map();
paramNames = keys(streamAnnos);

for i = 1:length(paramNames)
    key = paramNames{i};
    pidName = key;
    pid = streamAnnos(key);
    
    % Get the annotations associated with the pid
    if isnan(pid)
        pidAnnos = annotations(cellfun(@(x) isnan(str2double(x)), annotations.parameters), :);
    else
        pidAnnos = annotations(cellfun(@(x) str2double(x) == pid, annotations.parameters), :);
    end
    pidAnnos = sortrows(pidAnnos, 'qcFlag');
    
    % Create an array of flags (set to pass by default)
    pidFlags = ones(height(ds), 1);
    timeVec = ds.Properties.RowTimes;
    
    % For each annotation, set the qcFlag for each respective time period
    for ind = 1:height(pidAnnos)
        beginDT = pidAnnos.beginDate(ind);
        endDT = pidAnnos.endDate(ind);
        qcFlag = double(pidAnnos.qcFlag(ind));
        
        % Convert the time to actual datetimes
        beginDT = datetime(beginDT, "TimeZone", "UTC");
        if ismissing(endDT) || isempty(endDT{1})
            endDT = datetime('now', 'TimeZone', 'UTC');
        else
            endDT = datetime(endDT, "TimeZone", "UTC");
        end
        
        % Set the qcFlags for the given time range
        timeMask = (timeVec >= beginDT) & (timeVec <= endDT);
        pidFlags(timeMask) = qcFlag;
    end
    
    % Save the results
    flagsDict(pidName) = pidFlags;
end

% Add the flag results to the dataset
flagKeys = keys(flagsDict);
for i = 1:length(flagKeys)
    key = flagKeys{i};
    
    % Generate a variable name
    varName = [lower(key), '_annotations_qc_results'];
    varName = strrep(varName, '-', '_');  % Replace any hyphens with underscores
    
    % Build the comment
    if strcmpi(key, 'rollup')
        comment = ['These QC flags are a rollup summary which represents a ' ...
                  'Human-in-the-loop (HITL) assessment of the data quality ' ...
                  'for all variables in the dataset.'];
    else
        comment = sprintf(['These QC flags represent a Human-in-the-loop (HITL) ' ...
                         'assessment of the data quality for the specific data ' ...
                         'variable %s.'], key);
    end
        
    % Add to the timetable with custom properties for metadata
    ds.(varName) = flagsDict(key);
    ds.Properties.VariableDescriptions{varName} = comment;
    ds.Properties.VariableUnits{varName} = '';
end
