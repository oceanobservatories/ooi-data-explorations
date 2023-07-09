function fh = common()
% Function handles to allow for calling of the common functions across multiple
% processing functions
ooinet_defaults
fh.list_sites = @list_sites;
fh.list_nodes = @list_nodes;
fh.list_sensors = @list_sensors;
fh.list_methods = @list_methods;
fh.list_streams = @list_streams;
fh.list_deployments = @list_deployments;
end %function

% site, node and sensor information
function sites = list_sites()
% LIST_SITES Return a list of all available sites in the system.
%
%   Returns a list of all the available sites in the system. The list can then
%   be used to either iterate over the sites programmatically or inform the user
%   of the available sites and their codes.
%
% C. Wingard, 2020-01-13

% request a list of the nodes that are available
sites = webread([BASE_URL SENSOR_URL], options);
end %function

function nodes = list_nodes(site)
% LIST_NODES Based on the site name, list the nodes that are available.
%
%   Uses the site designator to create a list of the nodes that are available.
%   The returned list can then be used to either iterate over the nodes
%   programmatically or inform the user of the available nodes and their codes.
%
% C. Wingard, 2019-12-11

% request a list of the nodes that are available
nodes = webread([BASE_URL DEPLOY_URL upper(site)], options);
end %function

function sensors = list_sensors(site, node)
% LIST_SENSORS Based on the site and node name, list the sensors that are available.
%
%   Uses the site and node designators to create a list of the sensors that are
%   available. The returned list can then be used to either iterate over the
%   sensors programmatically or inform the user of the available sensors and
%   their codes.
%
% C. Wingard, 2019-12-11

% request a list of the nodes that are available
sensors = webread([BASE_URL DEPLOY_URL upper(site) '/' upper(node)], options);
end %function

function methods = list_methods(site, node, sensor)
% LIST_METHODS Based on the site, node and sensor name, list the data delivery methods that are available.
%
%   Uses the sitenode designators to create a list of the sensors that are
%   available. The returned list can then be used to either iterate over the
%   sensors programmatically or inform the user of the available sensors and
%   their codes.
%
% C. Wingard, 2019-12-11

% request a list of the nodes that are available
methods = webread([BASE_URL SENSOR_URL upper(site) '/' upper(node) '/' upper(sensor)], options);
end %function

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

% request a list of the nodes that are available
streams = webread([BASE_URL SENSOR_URL upper(site) '/' upper(node) '/' upper(sensor) '/' (method)], options);
end %function

function deployments = list_deployments(site, node, sensor)
% LIST_DEPLOYMENTS Based on the site, node and sensor name, list the available
% mooring deployment numbers for this sensor.
%
%   Uses the site, node and sensor designators to create a list of the mooring
%   deployments that are available for this specific sensor (combination of the
%   site, node and sensor). The returned list can then be used to either iterate
%   over the sensors programmatically or inform the user of the available
%   sensors and their codes.
%
% C. Wingard, 2019-12-11

% request a list of the nodes that are available
deployments = webread([BASE_URL DEPLOY_URL upper(site) '/' upper(node) '/' upper(sensor)], options);
end %function
