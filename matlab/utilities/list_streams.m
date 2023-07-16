function streams = list_streams(site, node, sensor, method)
% LIST_STREAMS Based on the site, node and sensor name and the data delivery
% method, list the data streams that are available.
%
%   Uses the sitenode designators to create a list of the sensors that are
%   available. The returned list can then be used to either iterate over the
%   sensors programmatically or inform the user of the available sensors and
%   their codes.
%
% C. Wingard, 2019-12-11

% load the default names and access credentials
ooinet_defaults

% request a list of the nodes that are available
streams = webread([BASE_URL SENSOR_URL upper(site) '/' upper(node) '/' upper(sensor) '/' (method)], options);

end %function
