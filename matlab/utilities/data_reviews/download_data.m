function download_data(pathname, site, node, sensor, method, stream, strt, stop)
% download_data.m
%
% Uses the python request_data.py module from the ooicgsn_data_tools repo to
% request and download the NetCDF files via the OOI M2M system based on the
% site, node, sensor, method, stream and starting and ending dates.
% 
% C. Wingard, 2021-08-18

% convert the start and end dates from Matlab datenums to the string format
% expected by the M2M system.
strt = datestr(strt, 'yyyy-mm-ddTHH:MM:SS.000Z');
stop = datestr(stop, 'yyyy-mm-ddTHH:MM:SS.000Z');

% setup the full, absolute file name and path
filename = [site '.' node '.' sensor '.' method '.' stream '.' strt '-' stop '.nc'];
fname = fullfile(pathname, filename);

% create the command string to run the python data request
cmdstr = ['python request_data.py -s ' upper(site) ' -n ' upper(node) ' -r ' upper(sensor) ...
    ' -m ' method ' -t ' stream ' -b ' strt ' -e ' stop ' -f ' fname];

% issue the data request
[status, result] = system(cmdstr);
if status == 0
    sprintf('data successfully requested and downloaded')
else
    sprintf('data request failed')
    sprintf(result)
end %if
