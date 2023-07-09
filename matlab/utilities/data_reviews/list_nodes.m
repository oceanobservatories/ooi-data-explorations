function nodes = list_nodes(site)
% LIST_NODES Based on the site name, list the nodes that are available.
%
%   Uses the site designator to create a list of the nodes that are available.
%   The returned list can then be used to either iterate over the nodes
%   programmatically or inform the user of the available nodes and their codes.
%
% C. Wingard, 2019-12-11

% load the default names and access credentials
ooinet_defaults

% request a list of the nodes that are available
nodes = webread([BASE_URL DEPLOY_URL upper(site)], options);

end %function
