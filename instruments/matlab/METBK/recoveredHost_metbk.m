%Written By Craig Risien on June 7, 2019 using Matlab2018a

%.. set login details
api_key = "OOIAPI-853A3LA6QI3L62";
api_token = "WYAN89W5X4Z0QZ";
%.. set mooring info and time period of interest
start_date='2018-01-01T00:00:00.000Z';
end_date='2018-12-31T23:59:59.000Z';
mooring_name = 'CE02SHSM';
node = 'BUOY';

%.. Explicitly construct UFrame dataset names
if strcmp(mooring_name,'CE02SHSM') && strcmp(node,'BUOY')
    uframe_dataset_name = 'CE02SHSM/SBD11/06-METBKA000/recovered_host/metbk_a_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'BUOY')
    uframe_dataset_name = 'CE07SHSM/SBD11/06-METBKA000/recovered_host/metbk_a_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'BUOY')
    uframe_dataset_name = 'CE04OSSM/SBD11/06-METBKA000/recovered_host/metbk_a_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'BUOY')
    uframe_dataset_name = 'CE09OSSM/SBD11/06-METBKA000/recovered_host/metbk_a_dcl_instrument_recovered';
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
strings_to_match = sprintf('.*/deployment.*%s.*\\.nc$', strrep(dataset_name, '/', '-'));
nc_urls = cellfun(@(x) regexp(x, strings_to_match, 'match'), nc_urls_all, 'UniformOutput', 0);
nc_urls(cellfun(@isempty, nc_urls)) = [];  % cell elements are themselves cells
nc_urls = string(nc_urls(:));  % string array

time_array=[];salinity_array=[];wtmp_array=[];uwnd_array=[];vwnd_array=[];bpress_array=[];
atmp_array=[];relhum_array=[];lrad_array=[];srad_array=[];precip_array=[];
heatflx_array=[];latnflx_array=[];netlirr_array=[];sensflx_array=[];uwater_array=[];vwater_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    time_array(length(time_array)+1:length(time_array)+length(data)) = data;clear data
    %Seawater Temperature (oC)
    data=ncread(char(nc_urls(i,:)),'sea_surface_temperature');
    wtmp_array(length(wtmp_array)+1:length(wtmp_array)+length(data)) = data;clear data
    %Seawater Salinity
    data=ncread(char(nc_urls(i,:)),'met_salsurf');
    salinity_array(length(salinity_array)+1:length(salinity_array)+length(data)) = data;clear data
    %Mean Wind Velocity - Eastward Relative to True North (m/s)
    data=ncread(char(nc_urls(i,:)),'met_windavg_mag_corr_east');
    uwnd_array(length(uwnd_array)+1:length(uwnd_array)+length(data)) = data;clear data
    %Mean Wind Velocity - Northward Relative to True North (m/s)
    data=ncread(char(nc_urls(i,:)),'met_windavg_mag_corr_north');
    vwnd_array(length(vwnd_array)+1:length(vwnd_array)+length(data)) = data;clear data
    %Barometric Pressure (mbar)
    data=ncread(char(nc_urls(i,:)),'barometric_pressure');
    bpress_array(length(bpress_array)+1:length(bpress_array)+length(data)) = data;clear data
    %Air Temperature (oC)
    data=ncread(char(nc_urls(i,:)),'air_temperature');
    atmp_array(length(atmp_array)+1:length(atmp_array)+length(data)) = data;clear data
    %Relative Humidity (%)
    data=ncread(char(nc_urls(i,:)),'relative_humidity');
    relhum_array(length(relhum_array)+1:length(relhum_array)+length(data)) = data;clear data
    %Longwave Irradiance (W m-2)
    data=ncread(char(nc_urls(i,:)),'longwave_irradiance');
    lrad_array(length(lrad_array)+1:length(lrad_array)+length(data)) = data;clear data
    %Shortwave Irradiance (W m-2)
    data=ncread(char(nc_urls(i,:)),'shortwave_irradiance');
    srad_array(length(srad_array)+1:length(srad_array)+length(data)) = data;clear data
    %Precipitation (mm)
    data=ncread(char(nc_urls(i,:)),'precipitation');
    precip_array(length(precip_array)+1:length(precip_array)+length(data)) = data;clear data
    %The total net upward heat flux data product HEATFLX_MINUTE_L2 (W m-2)
    data=ncread(char(nc_urls(i,:)),'met_heatflx_minute');
    heatflx_array(length(heatflx_array)+1:length(heatflx_array)+length(data)) = data;clear data
    %The upward latent heat flux data product LATNFLX_MINUTE_L2 (W m-2)
    data=ncread(char(nc_urls(i,:)),'met_latnflx_minute');
    latnflx_array(length(latnflx_array)+1:length(latnflx_array)+length(data)) = data;clear data
    %The net upward longwave irradiance NETLIRR_MINUTE_L2 (W m-2)
    data=ncread(char(nc_urls(i,:)),'met_netlirr_minute');
    netlirr_array(length(netlirr_array)+1:length(netlirr_array)+length(data)) = data;clear data
    %The net upward sensible heat flux SENSFLX_MINUTE_L2 (W m-2)
    data=ncread(char(nc_urls(i,:)),'met_sensflx_minute');
    sensflx_array(length(sensflx_array)+1:length(sensflx_array)+length(data)) = data;clear data
    %Eastward Mean Point Seawater Velocity (m/s)
    data=ncread(char(nc_urls(i,:)),'eastward_velocity');
    uwater_array(length(uwater_array)+1:length(uwater_array)+length(data)) = data;clear data
    %Northward Mean Point Seawater Velocity (m/s)
    data=ncread(char(nc_urls(i,:)),'northward_velocity');
    vwater_array(length(vwater_array)+1:length(vwater_array)+length(data)) = data;clear data
    
end

time_array=datenum(1900,1,1,0,0,0)+(time_array/60/60/24);

%Download NDBC Data (46050 and 46097, which is the WMO# for CE02SHSM)
%https://dods.ndbc.noaa.gov/thredds/catalog/data/stdmet/46050/catalog.html
time_46050=double(ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/46050/46050h2018.nc','time'));
time_46050=datenum(1970,1,1,0,0,0)+(time_46050/60/60/24);
sst_46050=squeeze(ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/46050/46050h2018.nc','sea_surface_temperature'));
wdir_46050=squeeze(ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/46050/46050h2018.nc','wind_dir'));
wspd_46050=squeeze(ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/46050/46050h2018.nc','wind_spd'));
%convert from spd and dir to u and v wnd
wdir_46050 = 90 - wdir_46050 + 180;
ind = wdir_46050 > 180;
wdir_46050(ind) = wdir_46050(ind) - 360;
ind = wdir_46050 < -180;
wdir_46050(ind) = wdir_46050(ind) + 360;
uwnd_46050 = wspd_46050.*cos(deg2rad(wdir_46050));
vwnd_46050 = wspd_46050.*sin(deg2rad(wdir_46050));

%https://dods.ndbc.noaa.gov/thredds/catalog/data/stdmet/46097/catalog.html
time_46097=double(ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/46097/46097h2018.nc','time'));
time_46097=datenum(1970,1,1,0,0,0)+(time_46097/60/60/24);
sst_46097=squeeze(ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/46097/46097h2018.nc','sea_surface_temperature'));
wdir_46097=squeeze(ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/46097/46097h2018.nc','wind_dir'));
wspd_46097=squeeze(ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/46097/46097h2018.nc','wind_spd'));
%convert from spd and dir to u and v wnd
wdir_46097 = 90 - wdir_46097 + 180;
ind = wdir_46097 > 180;
wdir_46097(ind) = wdir_46097(ind) - 360;
ind = wdir_46097 < -180;
wdir_46097(ind) = wdir_46097(ind) + 360;
uwnd_46097 = wspd_46097.*cos(deg2rad(wdir_46097));
vwnd_46097 = wspd_46097.*sin(deg2rad(wdir_46097));


%Download SCOW Climaotlogy Fields
vwnd_scow_regress_coefficients=ncread('https://wilson.coas.oregonstate.edu/thredds/dodsC/CIOSS/SCOW/Coeff/wind_meridional_coefficients.nc','regress_coefficients');
sst_scow_regress_coefficients=ncread('https://wilson.coas.oregonstate.edu/thredds/dodsC/CIOSS/SCOW/Coeff/avhrr_sst_coefficients.nc','regress_coefficients');
lat_scow=ncread('https://wilson.coas.oregonstate.edu/thredds/dodsC/CIOSS/SCOW/Coeff/wind_meridional_coefficients.nc','latitude');
lon_scow=ncread('https://wilson.coas.oregonstate.edu/thredds/dodsC/CIOSS/SCOW/Coeff/wind_meridional_coefficients.nc','longitude');
%Get the Regression Coefficients for the grid cell closest to CE02
ce02_lat=44.6393;
ce02_lon=-124.304+360;
[~,ind_lat] = min(abs(ce02_lat-lat_scow));
[~,ind_lon] = min(abs(ce02_lon-lon_scow));
vwnd_scow_regress_coefficients = squeeze(vwnd_scow_regress_coefficients(:,ind_lon,ind_lat));
sst_scow_regress_coefficients = squeeze(sst_scow_regress_coefficients(:,ind_lon,ind_lat));
%Calculate daily Climatologies 
new_t=0:1:364;
f=1/365;
vwnd_scow_seasonal_cycle=(vwnd_scow_regress_coefficients(1)+vwnd_scow_regress_coefficients(2)*sin(2*pi*f*new_t)+vwnd_scow_regress_coefficients(3)*cos(2*pi*f*new_t)+vwnd_scow_regress_coefficients(4)*sin(4*pi*f*new_t)+vwnd_scow_regress_coefficients(5)*cos(4*pi*f*new_t));
sst_scow_seasonal_cycle=(sst_scow_regress_coefficients(1)+sst_scow_regress_coefficients(2)*sin(2*pi*f*new_t)+sst_scow_regress_coefficients(3)*cos(2*pi*f*new_t)+sst_scow_regress_coefficients(4)*sin(4*pi*f*new_t)+vwnd_scow_regress_coefficients(5)*cos(4*pi*f*new_t));
scow_time=datenum(2018,1,1,12,0,0):datenum(2018,12,31,12,0,0);

%Plot some data
ticksx=datenum(2018,1,1):datenum(2019,1,1);
doy=str2num(datestr(ticksx,7));
ind=find(doy==1);

subplot(211)
plot(time_array,vwnd_array,'.k')
hold on
plot(time_46050,vwnd_46050,'.b')
plot(time_46097,vwnd_46097,'.r')
plot(scow_time,vwnd_scow_seasonal_cycle,'.m')
legend('CE02SHSM','46050','46097','Climaotology','AutoUpdate','off')
plot([datenum(2018,1,1):datenum(2019,1,1)],zeros(366,1),'k')
axis([datenum(2018,1,1) datenum(2019,1,1) -20 20])
xticks(ticksx(ind))
xticklabels(datestr(ticksx(ind)))
ylabel('m/s')
title('Meridional Wind')
%SST
subplot(212)
plot(time_array,wtmp_array,'.k')
hold on
plot(time_46050,sst_46050,'.b')
plot(time_46097,sst_46097,'.r')
plot(scow_time,sst_scow_seasonal_cycle,'.m')
legend('CE02SHSM','46050','46097','Climaotology','AutoUpdate','off')
axis([datenum(2018,1,1) datenum(2019,1,1) 8 20])
xticks(ticksx(ind))
xticklabels(datestr(ticksx(ind)))
ylabel('^oC')
title('Seawater Temperature')
