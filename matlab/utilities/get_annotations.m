function annotations = get_annotations(site, node, sensor)
% GET_ANNOTATIONS Based on the site, node and sensor name download the
% annotation records associated with this sensor.
%
%   Uses the site, node and sensor designators to obtain the annotation records
%   for the instrument. The annotations represent the HITL QC efforts on the
%   part of the OOI data teams, and as such provide a great deal of valuable
%   information about the instrument of interest.
%
% C. Wingard, 2020-02-12

% load the default names and access credentials
ooinet_defaults

% request a list of the nodes that are available
annos = webread([BASE_URL ANNO_URL 'find?beginDT=0&refdes=' upper(site) '-' upper(node) '-' upper(sensor)], options);
varNames = {'id', 'subsite', 'node', 'sensor', 'stream', 'method', ...
    'parameters', 'beginDate', 'endDate', 'exclusionFlag', 'qcFlag', ...
    'source', 'annotation'};
varTypes = {'double', 'string', 'string', 'string', 'string', 'string', ...
    'string', 'string', 'string', 'string', 'string', 'string', 'string'};
annotations = table('Size', [length(annos), 13], 'VariableTypes', varTypes, ...
    'VariableNames', varNames);
for i = 1:length(annos)
    annotations.id(i) = annos(i).id;
    annotations.subsite(i) = annos(i).subsite;
    if ~isempty(annos(i).node)
        annotations.node(i) = annos(i).node;
    end %if
    if ~isempty(annos(i).sensor)
        annotations.sensor(i) = annos(i).sensor;
    end %if
    if ~isempty(annos(i).stream)
        annotations.stream(i) = annos(i).stream;
    end %if
    if ~isempty(annos(i).method)
        annotations.method(i) = annos(i).method;
    end %if
    if ~isempty(annos(i).parameters)
        parms = sprintf('[%s]', join(string(annos(i).parameters), ','));
        annotations.parameters{i} = parms;
    end %if
    if ~isempty(annos(i).qcFlag)
        annotations.qcFlag(i) = annos(i).qcFlag;
    end %if
    beginDate = datenum([1970 1 1 0 0 (annos(i).beginDT / 1000)]);
    annotations.beginDate(i) = datestr(beginDate, 'yyyy-mm-ddTHH:MM:SS');
    if ~isempty(annos(i).endDT)
        endDate = datenum([1970 1 1 0 0 (annos(i).endDT / 1000)]);
        annotations.endDate(i) = datestr(endDate, 'yyyy-mm-ddTHH:MM:SS');
    end %if
    annotations.exclusionFlag(i) = annos(i).exclusionFlag;
    annotations.source(i) = annos(i).source;
    annotations.annotation(i) = annos(i).annotation;
end %for
clear annos varNames varTypes i beginDate endDate

end %function
