function info = get_parameter_information(parameter_id)
% GET_PARAMETER_INFORMATION Use the Parameter ID# to retrieve information
% about the parameter: units, sources, data product ID, comments, etc
%
% C. Wingard, 2026-01-06

% load the default names and access credentials
ooinet_defaults

% get the parameter information
info = webread([BASE_URL PARAMETER_URL parameter_id], options);
