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

% load the default names and access credentials
ooinet_defaults

% request a list of the nodes that are available
deployments = webread([BASE_URL DEPLOY_URL upper(site) '/' upper(node) '/' upper(sensor)], options);

end %function
