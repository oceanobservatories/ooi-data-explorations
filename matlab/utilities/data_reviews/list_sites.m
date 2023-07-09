function sites = list_sites()
% LIST_SITES Return a list of all available sites in the system.
%
%   Returns a of all the available sites in the system. The list can then be
%   used to either iterate over the sites programmatically or inform the user of
%   the available sites and their codes.
%
% C. Wingard, 2020-01-13

% load the default names and access credentials
ooinet_defaults

% request a list of the nodes that are available
sites = webread([BASE_URL SENSOR_URL], options);

end %function
