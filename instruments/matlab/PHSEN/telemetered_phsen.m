
%Written By Craig Risien on June 26, 2019 using Matlab2018a

%.. set login details
api_key = "OOIAPI-853A3LA6QI3L62";
api_token = "WYAN89W5X4Z0QZ";
%.. set mooring info and time period of interest
start_date='2017-01-01T00:00:00.000Z';
end_date='2018-12-31T23:59:59.000Z';
mooring_name = 'CE07SHSM';
node = 'NSIF'; %BUOY, or NSIF

%.. Explicitly construct UFrame dataset names
if strcmp(mooring_name,'CE01ISSM') && strcmp(node,'NSIF')
    uframe_dataset_name = 'CE01ISSM/RID16/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'NSIF')
    uframe_dataset_name = 'CE02SHSM/RID26/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'NSIF')
    uframe_dataset_name = 'CE04OSSM/RID26/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'NSIF')
    uframe_dataset_name = 'CE06ISSM/RID16/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'NSIF')
    uframe_dataset_name = 'CE07SHSM/RID26/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'NSIF')
    uframe_dataset_name = 'CE09OSSM/RID26/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'MFN')
    uframe_dataset_name = 'CE01ISSM/MFD35/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'MFN')
    uframe_dataset_name = 'CE06ISSM/MFD35/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'MFN')
    uframe_dataset_name = 'CE07SHSM/MFD35/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'MFN')
    uframe_dataset_name = 'CE09OSSM/MFD35/06-PHSEND000/telemetered/phsen_abcdef_dcl_instrument';
else
    error('Illegal mooring_name or node or combination thereof.');
end

%.. construct the data URL
m2m_url = "https://ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv/";
data_url = strcat(m2m_url, uframe_dataset_name);

options = weboptions("Username", api_key, "Password", api_token, "Timeout", 120);
m2m_response = webread(data_url, "beginDT", start_date, "endDT", end_date, options);

for ii = 1:1800
    try
        check_complete = webread(join([m2m_response.allURLs{2}, "status.txt"], "/"));
        if strip(check_complete) == "complete"
            fprintf("request completed\n");
            response_status = 'request complete';
            break
        else
            pause(1);
            response_status = 'request incomplete';
            if mod(ii,10)==0, fprintf('%u ', ii); end
        end
    catch
        pause(1);
        response_status = 'request incomplete';
        if mod(ii,10)==0, fprintf('%u ', ii); end
    end
end

% now put together the list of available NetCDF files from the THREDDS server
thredds_url = "https://opendap.oceanobservatories.org/thredds/dodsC";
catalog = webread(m2m_response.allURLs{1});
nclist = regexp(catalog, '<a href=''([^>]+.nc)''>', 'tokens'); % cell elements are themselves cells

if isempty(nclist)
    response_status = 'm2m request completed: no uframe netcdf files';
    disp(['uframe_m2m_status: ' response_status]);
    return
else
    response_status = 'successful';
    disp(['uframe_m2m_status: ' response_status]);
end

disp(['THREDDS Catalog URL: ' m2m_response.allURLs{1}])

nclist = [nclist{:}]';  % nclist is now a cell array of character vectors
nclist = cellfun(@(x) x(22:end), nclist, 'UniformOutput', 0); % drop 'catalog.html?dataset=' from URL
%.. create a URL mapping of the files
nc_urls_all = strcat(thredds_url, "/", nclist);  % string array

%.. prune the list of mapped files to our instrument of interest by eliminating
%.. files associated with other instruments used to calculate data products
strings_to_match = sprintf('.*/deployment.*%s.*\\.nc$', strrep(uframe_dataset_name, '/', '-'));
nc_urls = cellfun(@(x) regexp(x, strings_to_match, 'match'), nc_urls_all, 'UniformOutput', 0);
nc_urls(cellfun(@isempty, nc_urls)) = [];  % cell elements are themselves cells
nc_urls = string(nc_urls(:));  % string array

time_array=[];ph_array=[];temperature_array=[];salinity_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    time_array(length(time_array)+1:length(time_array)+length(data)) = data;clear data
    %PH
    data=ncread(char(nc_urls(i,:)),'phsen_abcdef_ph_seawater');
    ph_array(length(ph_array)+1:length(ph_array)+length(data)) = data;clear data
    %Temperature
    data=ncread(char(nc_urls(i,:)),'phsen_thermistor_temperature');
    temperature_array(length(temperature_array)+1:length(temperature_array)+length(data)) = data;clear data
    %Salinity
    data=ncread(char(nc_urls(i,:)),'practical_salinity');
    salinity_array(length(salinity_array)+1:length(salinity_array)+length(data)) = data;clear data
    
end

time_array=datenum(1900,1,1,0,0,0)+(time_array/60/60/24);

%Plot some data
ticksx=datenum(2018,1,1):datenum(2019,1,1);
doy=str2num(datestr(ticksx,7));
ind=find(doy==1);

subplot(311)
plot(time_array,ph_array,'.k')
axis([datenum(2018,1,1) datenum(2019,1,1) 7.5 9])
xticks(ticksx(ind))
xticklabels(datestr(ticksx(ind)))
title(strcat(mooring_name,{' '},node,{' '},'pH'))

subplot(312)
plot(time_array,temperature_array,'.k')
axis([datenum(2018,1,1) datenum(2019,1,1) 6 16])
xticks(ticksx(ind))
xticklabels(datestr(ticksx(ind)))
title(strcat(mooring_name,{' '},node,{' '},'Temperature'))
ylabel('^oC')

subplot(313)
plot(time_array,salinity_array,'.k')
axis([datenum(2018,1,1) datenum(2019,1,1) 26 34])
xticks(ticksx(ind))
xticklabels(datestr(ticksx(ind)))
title(strcat(mooring_name,{' '},node,{' '},'Salinity'))
