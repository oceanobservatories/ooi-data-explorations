function metadata = get_sensor_information(site, node, sensor, deploy)
% GET_SENSOR_INFORMATION Based on the site, node and sensor names and the
% deployment number, get the sensor metadata.
% 
%   Uses the metadata information available from the system for an instrument
%   deployment to obtain the asset and calibration information for the specified
%   sensor and deployment. This information is part of the instrument metadata
%   specific to that deployment.
%
% C. Wingard, 2020-01-13

% load the default names and access credentials
ooinet_defaults

% request the deployment information
metadata = webread([BASE_URL DEPLOY_URL upper(site) '/' upper(node) '/' upper(sensor) '/' sprintf('%d', deploy)], options);
metadata = metadata.sensor;

end %function
