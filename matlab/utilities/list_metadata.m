function metadata = list_metadata(site, node, sensor)
% LIST_METADATA Based on the site, node and sensor names, return a metadata dictionary
% with the times and parameters available for a sensor.
%
%   Uses the site, node and sensor names to create a metadata dictionary with
%   the parameters and times available for that sensor. The returned dictionary
%   can then be used to determine available time ranges (there may be gaps),
%   parameter names, units, variable types, etc.
%
% C. Wingard, 2020-01-13

% load the default names and access credentials
ooinet_defaults

% request a structured array with a subset of metadata available about the
% sensor of interest.
metadata = webread([BASE_URL SENSOR_URL upper(site) '/' upper(node) '/' upper(sensor) '/metadata'], options);

end %function
