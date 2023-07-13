function methods = list_methods(site, node, sensor)
% LIST_METHODS Based on the site, node and sensor name, list the data delivery methods that are available.
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
methods = webread([BASE_URL SENSOR_URL upper(site) '/' upper(node) '/' upper(sensor)], options);

end %function
