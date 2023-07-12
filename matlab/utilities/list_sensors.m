function sensors = list_sensors(site, node)
% LIST_SENSORS Based on the site and node name, list the sensors that are available.
%
%   Uses the site and node designators to create a list of the sensors that are
%   available. The returned list can then be used to either iterate over the
%   sensors programmatically or inform the user of the available sensors and
%   their codes.
%
% C. Wingard, 2019-12-11

% load the default names and access credentials
ooinet_defaults

% request a list of the nodes that are available
sensors = webread([BASE_URL DEPLOY_URL upper(site) '/' upper(node)], options);

end %function
