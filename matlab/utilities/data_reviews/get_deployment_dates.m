function [start, stop] = get_deployment_dates(site, node, sensor, deploy)
% GET_DEPLOYMENT_DATES Based on the site, node and sensor names and the 
% deployment number, determine the start and stop date/times for a deployment.
% 
%   Uses the metadata information available from the system for an instrument
%   deployment to determine the starting and ending dates for that deployment.
%   These dates (in Matlabs datenum format) can then be used to bound requests
%   to the system for specific date ranges specific to a deployment.
%
% C. Wingard, 2019-12-11

% load the default names and access credentials
ooinet_defaults

% request the deployment information
dates = webread([BASE_URL DEPLOY_URL upper(site) '/' upper(node) '/' upper(sensor) '/' sprintf('%d', deploy)], options);

% convert the timestamps from the unix epoch to a matlab datenum
if ~isempty(dates)
    start = datenum([1970 1 1 0 0 (dates.eventStartTime / 1000)]);
    if dates.eventStopTime
        stop = datenum([1970 1 1 0 0 (dates.eventStopTime / 1000)]);
    else
        stop = 0;
    end %if
else
    start = 0; stop = 0;
end %if

end %function
