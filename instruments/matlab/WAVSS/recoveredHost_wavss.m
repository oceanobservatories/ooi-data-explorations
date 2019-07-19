%%
%.. set login and URL details
api_key = "OOIAPI-853A3LA6QI3L62";
api_token = "WYAN89W5X4Z0QZ";
options = weboptions("Username", api_key, "Password", api_token, "Timeout", 120);

%.. set time period of interest
start_date='2015-01-01T00:00:00.000Z';
end_date='2019-4-30T23:59:59.000Z';

node = 'BUOY';
instrument_class = 'WAVSS';

%%
%Get CE04OSSM M2M URL
mooring_name = 'CE04OSSM';
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

ce04_time_array=[];ce04_wave_height=[];ce04_wave_period=[];ce04_wave_direction=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    ce04_time_array(length(ce04_time_array)+1:length(ce04_time_array)+length(data)) = data;clear data
    %Wave Direction
    data=ncread(char(nc_urls(i,:)),'mean_direction');
    ce04_wave_direction(length(ce04_wave_direction)+1:length(ce04_wave_direction)+length(data)) = data;clear data
    %Wave Period
    data=ncread(char(nc_urls(i,:)),'peak_wave_period');
    ce04_wave_period(length(ce04_wave_period)+1:length(ce04_wave_period)+length(data)) = data;clear data
    %Wave Height
    data=ncread(char(nc_urls(i,:)),'significant_wave_height');
    ce04_wave_height(length(ce04_wave_height)+1:length(ce04_wave_height)+length(data)) = data;clear data
    
end

ce04_time_array=datenum(1900,1,1,0,0,0)+(ce04_time_array/60/60/24);

%%
%Get CE02SHSM M2M URL
mooring_name = 'CE02SHSM';
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

ce02_time_array=[];ce02_wave_height=[];ce02_wave_period=[];ce02_wave_direction=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    ce02_time_array(length(ce02_time_array)+1:length(ce02_time_array)+length(data)) = data;clear data
    %Wave Direction
    data=ncread(char(nc_urls(i,:)),'mean_direction');
    ce02_wave_direction(length(ce02_wave_direction)+1:length(ce02_wave_direction)+length(data)) = data;clear data
    %Wave Period
    data=ncread(char(nc_urls(i,:)),'peak_wave_period');
    ce02_wave_period(length(ce02_wave_period)+1:length(ce02_wave_period)+length(data)) = data;clear data
    %Wave Height
    data=ncread(char(nc_urls(i,:)),'significant_wave_height');
    ce02_wave_height(length(ce02_wave_height)+1:length(ce02_wave_height)+length(data)) = data;clear data
    
end

ce02_time_array=datenum(1900,1,1,0,0,0)+(ce02_time_array/60/60/24);

%%
%Get CE01ISSM M2M URL
mooring_name = 'CE01ISSM';
node = 'MFN';
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

ce01_time_array=[];ce01_wave_height=[];ce01_wave_period=[];ce01_wave_direction=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    ce01_time_array(length(ce01_time_array)+1:length(ce01_time_array)+length(data)) = data;clear data
    %Wave Direction
    data=ncread(char(nc_urls(i,:)),'peak_wave_direction');
    ce01_wave_direction(length(ce01_wave_direction)+1:length(ce01_wave_direction)+length(data)) = data;clear data
    %Wave Period
    data=ncread(char(nc_urls(i,:)),'peak_wave_period');
    ce01_wave_period(length(ce01_wave_period)+1:length(ce01_wave_period)+length(data)) = data;clear data
    %Wave Height
    data=ncread(char(nc_urls(i,:)),'significant_wave_height');
    ce01_wave_height(length(ce01_wave_height)+1:length(ce01_wave_height)+length(data)) = data;clear data
    
end

ce01_time_array=datenum(1900,1,1,0,0,0)+(ce01_time_array/60/60/24);

%%
%Plot some data
for year = 2015:2019
    
    ticksx=datenum(year,1,1):datenum(year+1,1,1);
    doy=str2num(datestr(ticksx,7));
    ind=find(doy==1);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    
    subplot(311)
    plot(ce01_time_array,ce01_wave_height,'xk')
    hold on
    plot(ce02_time_array,ce02_wave_height,'.b')
    plot(ce04_time_array,ce04_wave_height,'.r')
    legend('CE01','CE02','CE04','location','northwest')
    axis([datenum(year,1,1) datenum(year+1,1,1) 0 15])
    xticks(ticksx(ind))
    xticklabels(datestr(ticksx(ind)))
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('meters')
    title('Significant Wave Height','fontweight','normal');
    
    subplot(312)
    plot(ce01_time_array,ce01_wave_period,'xk')
    hold on
    plot(ce02_time_array,ce02_wave_period,'.b')
    plot(ce04_time_array,ce04_wave_period,'.r')
    axis([datenum(year,1,1) datenum(year+1,1,1) 0 20])
    xticks(ticksx(ind))
    xticklabels(datestr(ticksx(ind)))
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('seconds')
    title('Wave Period','fontweight','normal');
    
    subplot(313)
    plot(ce01_time_array,ce01_wave_direction,'xk')
    hold on
    plot(ce02_time_array,ce02_wave_direction,'.b')
    plot(ce04_time_array,ce04_wave_direction,'.r')
    axis([datenum(year,1,1) datenum(year+1,1,1) 0 400])
    xticks(ticksx(ind))
    xticklabels(datestr(ticksx(ind)))
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('degress')
    title('Wave Direction','fontweight','normal');
    
end
