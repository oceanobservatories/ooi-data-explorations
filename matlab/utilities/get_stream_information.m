function definition = get_stream_information(stream)
% GET_STREAM_INFORMATION Use the stream name to retrieve information about the
% stream contents: parameters, units, sources, etc. 
%
%   Detailed explanation goes here
%
% C. Wingard, 2019-12-11

% load the default names and access credentials
ooinet_defaults

% request a list of the nodes that are available
definition = webread([BASE_URL STREAM_URL stream], options);

end %function
