function vocab = get_vocabulary(site, node, sensor)
% GET_VOCABULARY Based on the site, node and sensor name download the vocabulary
% record defining this sensor.
%
%   Uses the site, node and sensor designators to obtain the vocabulary record
%   for this sensor. This is a brief description of the instrument (vendor,
%   model, deployment site and depth) that may be useful from a metadata
%   standpoint.
%
% C. Wingard, 2019-12-11

% load the default names and access credentials
ooinet_defaults

% request a list of the nodes that are available
vocab = webread([BASE_URL VOCAB_URL upper(site) '/' upper(node) '/' upper(sensor)], options);

end %function
