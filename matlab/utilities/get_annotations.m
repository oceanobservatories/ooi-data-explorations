function annotations = get_annotations(site, varargin)
% GET_ANNOTATIONS Based on the site, node and sensor names.
%
%   Uses the site, node and sensor names to obtain the annotation records
%   for the instrument. Any request based on the site only, the site +
%   node, or the site + node + sensor will work, allowing the user to
%   request annotations in granular fashion for the site, node, sensor of
%   interest.
%
%   The annotations represent the HITL QC efforts on the part of the OOI
%   data teams, and as such provide a great deal of valuable information
%   about the instrument of interest.
%
%   Examples:
%
%   % All annotations for the site CE01ISSP
%   annotations = get_annotations('CE01ISSP');  
%
%   % All annotations for the site CE01ISSP and the node SP001
%   annotations = get_annotations('CE01ISSP', 'SP001');  
%
%   % All annotations for the site CE01ISSP, the node 'SP001', and the
%   sensor '09-CTDPFJ000'.
%   annotations = get_annotations('CE01ISSP', 'SP001', '09-CTDPFJ');
%
%   Note, the order of the site, node and sensor is important, and reflects
%   how the reference designators are constructed. If your request fails,
%   make sure you are properly structuring the request (use examples above,
%   which should work).
%
% C. Wingard, 2020-02-12 -- Original code
% C. Wingard, 2024-08-26 -- Updated to allow for granular requests

% load the default names and access credentials
ooinet_defaults

% request the annotations based on either the site, the site and node, or
% the site, node and sensor
switch length(varargin)
    case 0
        annos = webread([BASE_URL ANNO_URL 'find?beginDT=0&refdes=' upper(site)], options);
    case 1
        node = varargin{1};
        annos = webread([BASE_URL ANNO_URL 'find?beginDT=0&refdes=' upper(site) '-' upper(node)], options);
    case 2
        node = varargin{1};
        sensor = varargin{2};
        annos = webread([BASE_URL ANNO_URL 'find?beginDT=0&refdes=' upper(site) '-' upper(node) '-' upper(sensor)], options);
    otherwise
        error 'Incorrect numberof input arguements';
end %switch

% set up the annotations table
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
