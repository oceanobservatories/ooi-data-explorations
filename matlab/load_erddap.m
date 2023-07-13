function data = load_erddap(dataset_id, variables, times)
% LOAD_ERRDAP Loading data from the OOI Data Explorer ERDDAP Server
%
% Simple function to download date from the OOI Data Explorer ERDDAP server.
% Users will need to know the Dataset ID, the variable(s) of interest, and a
% time range to bound the size of the request.
%
% INPUTS:
%
%   dataset_id -- the ERDDAP ID of the dataset of interest
%   variables -- A string of the variables to download formatted as a 
%       comma-separated list (e.g., "temperature,pressure")
%   times -- start and end dates in an array formatted 'yyyy-mm-ddTHH:MM:SSZ'
%
% OUTPUTS:
%
%   data -- data downloaded from the ERDDAP server formatted as a timetable
%
% C. Wingard, 2021-12-20

% setup the request URL to download data from the Data Explorer ERDDAP
erddap_url = 'http://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap';
url = [erddap_url '/' dataset_id '.mat?time,' variables '&time>=' char(times(1, :)) '&time<=' char(times(2, :))];

% download the data as a MAT file and load it
websave('temp.mat', url);
m = matfile('temp.mat');
data = m.(replace(dataset_id, '-', '_'));
delete('temp.mat')
clear url m

% convert the time record from seconds since 1970 to a Matlab datenum ...
%data.time = data.time / 86400 + 719529;

% ... or a datetime object
data.time = datetime(data.time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

% convert the structured array to a timetable
data = table2timetable(struct2table(data), "RowTimes", "time");
end %function
