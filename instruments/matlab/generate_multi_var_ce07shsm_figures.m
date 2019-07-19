%%
clear all
close all

%.. set login and URL details
api_key = "OOIAPI-853A3LA6QI3L62";
api_token = "WYAN89W5X4Z0QZ";
options = weboptions("Username", api_key, "Password", api_token, "Timeout", 120);

%.. set mooring and time period of interest
mooring_name = 'CE07SHSM';
start_date='2015-01-01T00:00:00.000Z';
end_date='2019-06-30T23:59:59.000Z';

%%
%Get METBK Data
node = 'BUOY';
instrument_class = 'METBK';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

metbk_time_array=[];metbk_salinity_array=[];metbk_water_temp_array=[];metbk_uwnd_array=[];metbk_vwnd_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    metbk_time_array(length(metbk_time_array)+1:length(metbk_time_array)+length(data)) = data;clear data
    %Seawater Temperature (oC)
    data=ncread(char(nc_urls(i,:)),'sea_surface_temperature');
    metbk_water_temp_array(length(metbk_water_temp_array)+1:length(metbk_water_temp_array)+length(data)) = data;clear data
    %Seawater Salinity
    data=ncread(char(nc_urls(i,:)),'met_salsurf');
    metbk_salinity_array(length(metbk_salinity_array)+1:length(metbk_salinity_array)+length(data)) = data;clear data
    %Mean Wind Velocity - Eastward Relative to True North (m/s)
    data=ncread(char(nc_urls(i,:)),'met_windavg_mag_corr_east');
    metbk_uwnd_array(length(metbk_uwnd_array)+1:length(metbk_uwnd_array)+length(data)) = data;clear data
    %Mean Wind Velocity - Northward Relative to True North (m/s)
    data=ncread(char(nc_urls(i,:)),'met_windavg_mag_corr_north');
    metbk_vwnd_array(length(metbk_vwnd_array)+1:length(metbk_vwnd_array)+length(data)) = data;clear data
    
end

metbk_time_array=datenum(1900,1,1,0,0,0)+(metbk_time_array/60/60/24);

%%
%Get NSIF CTD data
node = 'NSIF';
instrument_class = 'CTDBP';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

ctd_nsif_time_array=[];ctd_nsif_salinity_array=[];ctd_nsif_temperature_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    ctd_nsif_time_array(length(ctd_nsif_time_array)+1:length(ctd_nsif_time_array)+length(data)) = data;clear data
    %Seawater Temperature (oC)
    data=ncread(char(nc_urls(i,:)),'temp');
    ctd_nsif_temperature_array(length(ctd_nsif_temperature_array)+1:length(ctd_nsif_temperature_array)+length(data)) = data;clear data
    %Seawater Salinity
    data=ncread(char(nc_urls(i,:)),'practical_salinity');
    ctd_nsif_salinity_array(length(ctd_nsif_salinity_array)+1:length(ctd_nsif_salinity_array)+length(data)) = data;clear data
    
end

ctd_nsif_time_array=datenum(1900,1,1,0,0,0)+(ctd_nsif_time_array/60/60/24);

%%
%Get MFN CTD data
node = 'MFN';
instrument_class = 'CTDBP';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

ctd_mfn_time_array=[];ctd_mfn_salinity_array=[];ctd_mfn_temperature_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    ctd_mfn_time_array(length(ctd_mfn_time_array)+1:length(ctd_mfn_time_array)+length(data)) = data;clear data
    %Seawater Temperature (oC)
    data=ncread(char(nc_urls(i,:)),'temp');
    ctd_mfn_temperature_array(length(ctd_mfn_temperature_array)+1:length(ctd_mfn_temperature_array)+length(data)) = data;clear data
    %Seawater Salinity
    data=ncread(char(nc_urls(i,:)),'practical_salinity');
    ctd_mfn_salinity_array(length(ctd_mfn_salinity_array)+1:length(ctd_mfn_salinity_array)+length(data)) = data;clear data
    
end

ctd_mfn_time_array=datenum(1900,1,1,0,0,0)+(ctd_mfn_time_array/60/60/24);

%%
%Get PCO2A Data
node = 'BUOY';
instrument_class = 'PCO2A';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

pco2a_time_array=[];pco2a_buoy_array=[];pco2w_buoy_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    pco2a_time_array(length(pco2a_time_array)+1:length(pco2a_time_array)+length(data)) = data;clear data
    %PCO2A
    data=ncread(char(nc_urls(i,:)),'partial_pressure_co2_atm');
    pco2a_buoy_array(length(pco2a_buoy_array)+1:length(pco2a_buoy_array)+length(data)) = data;clear data
    %PCO2W
    data=ncread(char(nc_urls(i,:)),'partial_pressure_co2_ssw');
    pco2w_buoy_array(length(pco2w_buoy_array)+1:length(pco2w_buoy_array)+length(data)) = data;clear data
    
end

pco2a_time_array=datenum(1900,1,1,0,0,0)+(pco2a_time_array/60/60/24);

%%
%Get PCO2W Data
node = 'MFN';
instrument_class = 'PCO2W';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

pco2w_mfn_time_array=[];pco2w_mfn_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    pco2w_mfn_time_array(length(pco2w_mfn_time_array)+1:length(pco2w_mfn_time_array)+length(data)) = data;clear data
    %PCO2
    data=ncread(char(nc_urls(i,:)),'pco2_seawater');
    pco2w_mfn_array(length(pco2w_mfn_array)+1:length(pco2w_mfn_array)+length(data)) = data;clear data
       
end

pco2w_mfn_time_array=datenum(1900,1,1,0,0,0)+(pco2w_mfn_time_array/60/60/24);

%%
%Get NSIF PHSEN Data
node = 'NSIF';
instrument_class = 'PHSEN';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

phsen_nsif_time_array=[];phsen_nsif_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    phsen_nsif_time_array(length(phsen_nsif_time_array)+1:length(phsen_nsif_time_array)+length(data)) = data;clear data
    %PH
    data=ncread(char(nc_urls(i,:)),'phsen_abcdef_ph_seawater');
    phsen_nsif_array(length(phsen_nsif_array)+1:length(phsen_nsif_array)+length(data)) = data;clear data
    
end

phsen_nsif_time_array=datenum(1900,1,1,0,0,0)+(phsen_nsif_time_array/60/60/24);

%%
%Get MFN PHSEN Data
node = 'MFN';
instrument_class = 'PHSEN';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

phsen_mfn_time_array=[];phsen_mfn_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    phsen_mfn_time_array(length(phsen_mfn_time_array)+1:length(phsen_mfn_time_array)+length(data)) = data;clear data
    %PH
    data=ncread(char(nc_urls(i,:)),'phsen_abcdef_ph_seawater');
    phsen_mfn_array(length(phsen_mfn_array)+1:length(phsen_mfn_array)+length(data)) = data;clear data
    
end

phsen_mfn_time_array=datenum(1900,1,1,0,0,0)+(phsen_mfn_time_array/60/60/24);

%%
%Get NSIF DO Data
node = 'NSIF';
instrument_class = 'DOSTA';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

dosta_nsif_time_array=[];dosta_nsif_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    dosta_nsif_time_array(length(dosta_nsif_time_array)+1:length(dosta_nsif_time_array)+length(data)) = data;clear data
    %L2 Pressure Temp Sal Corrected
    data=ncread(char(nc_urls(i,:)),'dissolved_oxygen'); 
    dosta_nsif_array(length(dosta_nsif_array)+1:length(dosta_nsif_array)+length(data)) = data;clear data
    
end

dosta_nsif_time_array=datenum(1900,1,1,0,0,0)+(dosta_nsif_time_array/60/60/24);

%%
%Get MFN DO Data
node = 'MFN';
instrument_class = 'DOSTA';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

dosta_mfn_time_array=[];dosta_mfn_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    dosta_mfn_time_array(length(dosta_mfn_time_array)+1:length(dosta_mfn_time_array)+length(data)) = data;clear data
    %L2 Pressure Temp Sal Corrected
    data=ncread(char(nc_urls(i,:)),'dissolved_oxygen'); 
    dosta_mfn_array(length(dosta_mfn_array)+1:length(dosta_mfn_array)+length(data)) = data;clear data
    
end

dosta_mfn_time_array=datenum(1900,1,1,0,0,0)+(dosta_mfn_time_array/60/60/24);

%%
%Get NSIF FLORT Data
node = 'NSIF';
instrument_class = 'FLORT';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Call
[nc_urls] = M2M_Call(uframe_dataset_name,data_url,start_date,end_date,options);

flort_nsif_time_array=[];flort_nsif_chla_array=[];flort_nsif_cdom_array=[];flort_nsif_backscatter_array=[];

for i = 1:length(nc_urls)
    
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    flort_nsif_time_array(length(flort_nsif_time_array)+1:length(flort_nsif_time_array)+length(data)) = data;clear data
    %CHLA
    data=ncread(char(nc_urls(i,:)),'fluorometric_chlorophyll_a');
    flort_nsif_chla_array(length(flort_nsif_chla_array)+1:length(flort_nsif_chla_array)+length(data)) = data;clear data
    %CDOM
    data=ncread(char(nc_urls(i,:)),'fluorometric_cdom');
    flort_nsif_cdom_array(length(flort_nsif_cdom_array)+1:length(flort_nsif_cdom_array)+length(data)) = data;clear data
    %Backscatter
    data=ncread(char(nc_urls(i,:)),'optical_backscatter');
    flort_nsif_backscatter_array(length(flort_nsif_backscatter_array)+1:length(flort_nsif_backscatter_array)+length(data)) = data;clear data
    
end

flort_nsif_time_array=datenum(1900,1,1,0,0,0)+(flort_nsif_time_array/60/60/24);

%%
%Get NSIF NUTNR Data from all three recoveredHost streams
node = 'NSIF';
instrument_class = 'NUTNR';

%Get M2M URL
[uframe_dataset_name,data_url]=RecoveredHost_M2M_URLs(mooring_name,node,instrument_class);

%Make M2M Calls
nutnr_time_array=[];nutnr_nitrate_concentration_array=[];
[nc_urls] = M2M_Call(uframe_dataset_name{1,1},data_url{1,1},start_date,end_date,options);
for i = 1:length(nc_urls)
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    nutnr_time_array(length(nutnr_time_array)+1:length(nutnr_time_array)+length(data)) = data;clear data
    %NO3
    data=ncread(char(nc_urls(i,:)),'nitrate_concentration');
    nutnr_nitrate_concentration_array(length(nutnr_nitrate_concentration_array)+1:length(nutnr_nitrate_concentration_array)+length(data)) = data;clear data
end
[nc_urls] = M2M_Call(uframe_dataset_name{1,2},data_url{1,2},start_date,end_date,options);
for i = 1:length(nc_urls)
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    nutnr_time_array(length(nutnr_time_array)+1:length(nutnr_time_array)+length(data)) = data;clear data
    %NO3
    data=ncread(char(nc_urls(i,:)),'salinity_corrected_nitrate');
    nutnr_nitrate_concentration_array(length(nutnr_nitrate_concentration_array)+1:length(nutnr_nitrate_concentration_array)+length(data)) = data;clear data
end
[nc_urls] = M2M_Call(uframe_dataset_name{1,3},data_url{1,3},start_date,end_date,options);
for i = 1:length(nc_urls)
    %Time (seconds since 1900-01-01 0:0:0)
    data=ncread(char(nc_urls(i,:)),'time');
    nutnr_time_array(length(nutnr_time_array)+1:length(nutnr_time_array)+length(data)) = data;clear data
    %NO3
    data=ncread(char(nc_urls(i,:)),'salinity_corrected_nitrate');
    nutnr_nitrate_concentration_array(length(nutnr_nitrate_concentration_array)+1:length(nutnr_nitrate_concentration_array)+length(data)) = data;clear data
end

nutnr_time_array=datenum(1900,1,1,0,0,0)+(nutnr_time_array/60/60/24);

%%
%Bin data hourly average
bin_edges = (datenum(2015,3,30,23,30,0)*24*60:60:datenum(2019,6,1,0,30,0)*24*60)/24/60; % make sure bounds outside actual data ranges

%NSIF NUTNR
ind = discretize(nutnr_time_array,bin_edges);
nutnr_nitrate_concentration_binned = accumarray(ind',nutnr_nitrate_concentration_array,[],@nanmedian);
nutnr_nitrate_concentration_binned(nutnr_nitrate_concentration_binned==0)=nan;
%NSIF NUTNR Time
nutnr_time_binned = accumarray(ind',nutnr_time_array,[],@nanmedian);
nutnr_time_binned(isnan(nutnr_nitrate_concentration_binned))=nan;clear ind

%METBK VWnd
ind = discretize(metbk_time_array,bin_edges);
metbk_vwnd_binned = accumarray(ind',metbk_vwnd_array,[],@nanmedian);
metbk_vwnd_binned(metbk_vwnd_binned==0)=nan;
metbk_vwnd_binned(metbk_vwnd_binned<-50)=nan;
metbk_vwnd_binned(metbk_vwnd_binned>50)=nan;
%METBK SST
metbk_water_temp_binned = accumarray(ind',metbk_water_temp_array,[],@nanmedian);
metbk_water_temp_binned(metbk_water_temp_binned==0)=nan;
%METBK SSS
metbk_salinity_binned = accumarray(ind',metbk_salinity_array,[],@nanmedian);
metbk_salinity_binned(metbk_salinity_binned==0)=nan;
%METBK Time
metbk_time_binned = accumarray(ind',metbk_time_array,[],@nanmedian);
metbk_time_binned(isnan(metbk_vwnd_binned))=nan;clear ind

%NSIF CTDBP WTemp
ind = discretize(ctd_nsif_time_array,bin_edges);
ctd_nsif_temperature_binned = accumarray(ind',ctd_nsif_temperature_array,[],@nanmedian);
ctd_nsif_temperature_binned(ctd_nsif_temperature_binned==0)=nan;
%NSIF CTDBP Salinity
ctd_nsif_salinity_binned = accumarray(ind',ctd_nsif_salinity_array,[],@nanmedian);
ctd_nsif_salinity_binned(ctd_nsif_salinity_binned==0)=nan;
%NSIF CTDBP Time
ctd_nsif_time_binned = accumarray(ind',ctd_nsif_time_array,[],@nanmedian);
ctd_nsif_time_binned(isnan(ctd_nsif_temperature_binned))=nan;clear ind

%MFN CTDBP WTemp
ind = discretize(ctd_mfn_time_array,bin_edges);
ctd_mfn_temperature_binned = accumarray(ind',ctd_mfn_temperature_array,[],@nanmedian);
ctd_mfn_temperature_binned(ctd_mfn_temperature_binned==0)=nan;
%MFN CTDBP Salinity
ctd_mfn_salinity_binned = accumarray(ind',ctd_mfn_salinity_array,[],@nanmedian);
ctd_mfn_salinity_binned(ctd_mfn_salinity_binned==0)=nan;
%MFN CTDBP Time
ctd_mfn_time_binned = accumarray(ind',ctd_mfn_time_array,[],@nanmedian);
ctd_mfn_time_binned(isnan(ctd_mfn_temperature_binned))=nan;clear ind

%BUOY PCO2A
ind = discretize(pco2a_time_array,bin_edges);
pco2a_buoy_binned = accumarray(ind',pco2a_buoy_array,[],@nanmedian);
pco2a_buoy_binned(pco2a_buoy_binned==0)=nan;
%BUOY PCO2W
pco2w_buoy_binned = accumarray(ind',pco2w_buoy_array,[],@nanmedian);
pco2w_buoy_binned(pco2w_buoy_binned==0)=nan;
%BUOY PCO2W Time
pco2a_time_binned = accumarray(ind',pco2a_time_array,[],@nanmedian);
pco2a_time_binned(isnan(pco2w_buoy_binned))=nan;clear ind

%MFN PCO2W
ind = discretize(pco2w_mfn_time_array,bin_edges);
pco2w_mfn_binned = accumarray(ind',pco2w_mfn_array,[],@nanmedian);
pco2w_mfn_binned(pco2w_mfn_binned==0)=nan;
%MFN PCO2W Time
pco2w_mfn_time_binned = accumarray(ind',pco2w_mfn_time_array,[],@nanmedian);
pco2w_mfn_time_binned(isnan(pco2w_mfn_binned))=nan;clear ind

%NSIF PHSEN
ind = discretize(phsen_nsif_time_array,bin_edges);
phsen_nsif_binned = accumarray(ind',phsen_nsif_array,[],@nanmedian);
phsen_nsif_binned(phsen_nsif_binned==0)=nan;
%NSIF PHSEN Time
phsen_nsif_time_binned = accumarray(ind',phsen_nsif_time_array,[],@nanmedian);
phsen_nsif_time_binned(isnan(phsen_nsif_binned))=nan;clear ind

%MFN PHSEN
ind = discretize(phsen_mfn_time_array,bin_edges);
phsen_mfn_binned = accumarray(ind',phsen_mfn_array,[],@nanmedian);
phsen_mfn_binned(phsen_mfn_binned==0)=nan;
%MFN PHSEN Time
phsen_mfn_time_binned = accumarray(ind',phsen_mfn_time_array,[],@nanmedian);
phsen_mfn_time_binned(isnan(phsen_mfn_binned))=nan;clear ind

%NSIF DOSTA
ind = discretize(dosta_nsif_time_array,bin_edges);
dosta_nsif_binned = accumarray(ind',dosta_nsif_array,[],@nanmedian);
dosta_nsif_binned(dosta_nsif_binned==0)=nan;
%NSIF DOSTA Time
dosta_nsif_time_binned = accumarray(ind',dosta_nsif_time_array,[],@nanmedian);
dosta_nsif_time_binned(isnan(dosta_nsif_binned))=nan;clear ind

%MFN DOSTA
ind = discretize(dosta_mfn_time_array,bin_edges);
dosta_mfn_binned = accumarray(ind',dosta_mfn_array,[],@nanmedian);
dosta_mfn_binned(dosta_mfn_binned==0)=nan;
%MFN DOSTA Time
dosta_mfn_time_binned = accumarray(ind',dosta_mfn_time_array,[],@nanmedian);
dosta_mfn_time_binned(isnan(dosta_mfn_binned))=nan;clear ind

%NSIF FLORT CHLA
ind = discretize(flort_nsif_time_array,bin_edges);
flort_nsif_chla_binned = accumarray(ind',flort_nsif_chla_array,[],@nanmedian);
flort_nsif_chla_binned(flort_nsif_chla_binned==0)=nan;
%NSIF FLORT CDOM
flort_nsif_cdom_binned = accumarray(ind',flort_nsif_cdom_array,[],@nanmedian);
flort_nsif_cdom_binned(flort_nsif_cdom_binned==0)=nan;
%NSIF FLORT BBACK
flort_nsif_backscatter_binned = accumarray(ind',flort_nsif_backscatter_array,[],@nanmedian);
flort_nsif_backscatter_binned(flort_nsif_backscatter_binned==0)=nan;
%NSIF FLORT Time
flort_nsif_time_binned = accumarray(ind',flort_nsif_time_array,[],@nanmedian);
flort_nsif_time_binned(isnan(flort_nsif_chla_binned))=nan;clear ind

%%
%Plot data March - October
for year = 2015:2018
    
    ticksx=datenum(year,3,1):datenum(year,11,1);
    doy=str2num(datestr(ticksx,7));
    ind=find(doy==1);
    
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    
    subplot('Position',[0.1 0.84 0.8 0.13])
    plot(metbk_time_binned,metbk_vwnd_binned,'k','linewidth',.8)
    hold on
    plot(datenum(year,3,1):datenum(year,11,1),zeros(length(datenum(year,3,1):datenum(year,11,1)),1),'k','linewidth',.4)
    axis([datenum(year,3,1) datenum(year,11,1) -15 15])
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('m/s')
    title('Meridional Winds','fontweight','normal','position',[datenum(year,3,11), 10.5, 0]);
    text(datenum(year,6,24,12,0,0),17,mooring_name,'fontsize',11)
    
    subplot('Position',[0.1 0.71 0.8 0.13])
    yyaxis left
    plot(metbk_time_binned,metbk_water_temp_binned,'-c','linewidth',.8)
    hold on
    plot(ctd_nsif_time_binned,ctd_nsif_temperature_binned,'-b')
    axis([datenum(year,3,1) datenum(year,11,1) 5 20])
    ylabel('DegC')
    ax = gca;ax.YColor = 'b';
    yyaxis right
    h=plot(metbk_time_binned,metbk_salinity_binned,'-k','linewidth',.8);set(h,'color',[1 .55 0])
    plot(ctd_nsif_time_binned,ctd_nsif_salinity_binned,'-r')
    axis([datenum(year,3,1) datenum(year,11,1) 20 35])
    ylabel('PSU')
    ax = gca;ax.YColor = 'r';
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    title('BUOY & NSIF Temperature and Salinity','fontweight','normal','position',[datenum(year,3,21), 18, 0]);
    legend('BUOY Tmp','NSIF Tmp','BUOY Sal','NSIF Sal')
    
    subplot('Position',[0.1 0.58 0.8 0.13])
    yyaxis left
    plot(pco2a_time_binned,pco2a_buoy_binned,'-k','linewidth',.8)
    axis([datenum(year,3,1) datenum(year,11,1) 0 500])
    ylabel('PCO2A (uatm)')
    ax = gca;ax.YColor = 'k';
    hold on
    yyaxis right
    h=plot(pco2a_time_binned,pco2w_buoy_binned,'-b','linewidth',.8);
    axis([datenum(year,3,1) datenum(year,11,1) 0 500])
    ylabel('PCO2W (uatm)')
    ax = gca;ax.YColor = 'b';
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    title('Surface PCO2','fontweight','normal','position',[datenum(year,3,10), 420, 0]);
    
    subplot('Position',[0.1 0.45 0.8 0.13])
    plot(phsen_nsif_time_binned,phsen_nsif_binned,'k','linewidth',.8)
    axis([datenum(year,3,1) datenum(year,11,1) 6 10])
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('pH')
    title('NSIF pH','fontweight','normal','position',[datenum(year,3,7), 9.4, 0]);
    
    subplot('Position',[0.1 0.32 0.8 0.13])
    plot(dosta_nsif_time_binned,dosta_nsif_binned,'k','linewidth',.8)
    axis([datenum(year,3,1) datenum(year,11,1) 0 600])
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('umol/kg')
    title('NSIF DO','fontweight','normal','position',[datenum(year,3,7), 510, 0]);
    
    subplot('Position',[0.1 0.19 0.8 0.13])
    yyaxis left
    plot(flort_nsif_time_binned,flort_nsif_chla_binned,'-k','linewidth',.8)
    axis([datenum(year,3,1) datenum(year,11,1) 0 30])
    ylabel('Chla (ug/L)')
    ax = gca;ax.YColor = 'k';
    hold on
    yyaxis right
    h=plot(flort_nsif_time_binned,flort_nsif_backscatter_binned,'-b','linewidth',.8);
    axis([datenum(year,3,1) datenum(year,11,1) 0 .1])
    ylabel('Backscatter (m-1)')
    ax = gca;ax.YColor = 'b';
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    title('NSIF Fluorometer','fontweight','normal','position',[datenum(year,3,10), 25, 0]);
    
    subplot('Position',[0.1 0.06 0.8 0.13])
    plot(nutnr_time_binned,nutnr_nitrate_concentration_binned,'k','linewidth',.8)
    hold on
    plot(datenum(year,3,1):datenum(year,11,1),zeros(length(datenum(year,3,1):datenum(year,11,1)),1),'k','linewidth',.4)
    axis([datenum(year,3,1) datenum(year,11,1) -5 30])
    xticks(ticksx(ind))
    xticklabels(datestr(ticksx(ind)))
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('umol/L')
    title('NSIF Nitrate','fontweight','normal','position',[datenum(year,3,7), 24, 0]);
    
    %save figure
    set(gcf,'color','w');
    export_fig('-r300',strcat(mooring_name,'_March_October_',num2str(year)),'-pdf')
end

%%
%Plot data September - April
for year = 2015:2018
    
    ticksx=datenum(year,9,1):datenum(year+1,5,1);
    doy=str2num(datestr(ticksx,7));
    ind=find(doy==1);
    
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    
    subplot('Position',[0.1 0.84 0.8 0.13])
    plot(metbk_time_binned,metbk_vwnd_binned,'k','linewidth',.8)
    hold on
    plot(datenum(year,9,1):datenum(year+1,5,1),zeros(length(datenum(year,9,1):datenum(year+1,5,1)),1),'k','linewidth',.4)
    axis([datenum(year,9,1) datenum(year+1,5,1) -15 15])
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('m/s')
    title('Meridional Winds','fontweight','normal','position',[datenum(year,9,11), 10.5, 0]);
    text(datenum(year,12,24,12,0,0),17,mooring_name,'fontsize',11)
    
    subplot('Position',[0.1 0.71 0.8 0.13])
    yyaxis left
    plot(metbk_time_binned,metbk_water_temp_binned,'-c','linewidth',.8)
    hold on
    plot(ctd_nsif_time_binned,ctd_nsif_temperature_binned,'-b')
    axis([datenum(year,9,1) datenum(year+1,5,1) 5 20])
    ylabel('DegC')
    ax = gca;ax.YColor = 'b';
    yyaxis right
    h=plot(metbk_time_binned,metbk_salinity_binned,'-k','linewidth',.8);set(h,'color',[1 .55 0])
    plot(ctd_nsif_time_binned,ctd_nsif_salinity_binned,'-r')
    axis([datenum(year,9,1) datenum(year+1,5,1) 20 35])
    ylabel('PSU')
    ax = gca;ax.YColor = 'r';
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    title('BUOY & NSIF Temperature and Salinity','fontweight','normal','position',[datenum(year,9,21), 18, 0]);
    legend('BUOY Tmp','NSIF Tmp','BUOY Sal','NSIF Sal')
    
    subplot('Position',[0.1 0.58 0.8 0.13])
    yyaxis left
    plot(pco2a_time_binned,pco2a_buoy_binned,'-k','linewidth',.8)
    axis([datenum(year,9,1) datenum(year+1,5,1) 0 500])
    ylabel('PCO2A (uatm)')
    ax = gca;ax.YColor = 'k';
    hold on
    yyaxis right
    h=plot(pco2a_time_binned,pco2w_buoy_binned,'-b','linewidth',.8);
    axis([datenum(year,9,1) datenum(year+1,5,1) 0 500])
    ylabel('PCO2W (uatm)')
    ax = gca;ax.YColor = 'b';
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    title('Surface PCO2','fontweight','normal','position',[datenum(year,9,10), 420, 0]);
    
    subplot('Position',[0.1 0.45 0.8 0.13])
    plot(phsen_nsif_time_binned,phsen_nsif_binned,'k','linewidth',.8)
    axis([datenum(year,9,1) datenum(year+1,5,1) 6 10])
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('pH')
    title('NSIF pH','fontweight','normal','position',[datenum(year,9,7), 9.4, 0]);
    
    subplot('Position',[0.1 0.32 0.8 0.13])
    plot(dosta_nsif_time_binned,dosta_nsif_binned,'k','linewidth',.8)
    axis([datenum(year,9,1) datenum(year+1,5,1) 0 600])
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('umol/kg')
    title('NSIF DO','fontweight','normal','position',[datenum(year,9,7), 510, 0]);
    
    subplot('Position',[0.1 0.19 0.8 0.13])
    yyaxis left
    plot(flort_nsif_time_binned,flort_nsif_chla_binned,'-k','linewidth',.8)
    axis([datenum(year,9,1) datenum(year+1,5,1) 0 30])
    ylabel('Chla (ug/L)')
    ax = gca;ax.YColor = 'k';
    hold on
    yyaxis right
    h=plot(flort_nsif_time_binned,flort_nsif_backscatter_binned,'-b','linewidth',.8);
    axis([datenum(year,9,1) datenum(year+1,5,1) 0 .1])
    ylabel('Backscatter (m-1)')
    ax = gca;ax.YColor = 'b';
    xticks(ticksx(ind))
    xticklabels([])
    set(gca,'TickLength',[0.005, 0.005])
    title('NSIF Fluorometer','fontweight','normal','position',[datenum(year,9,10), 25, 0]);
    
    subplot('Position',[0.1 0.06 0.8 0.13])
    plot(nutnr_time_binned,nutnr_nitrate_concentration_binned,'k','linewidth',.8)
    hold on
    plot(datenum(year,9,1):datenum(year+1,5,1),zeros(length(datenum(year,9,1):datenum(year+1,5,1)),1),'k','linewidth',.4)
    axis([datenum(year,9,1) datenum(year+1,5,1) -5 30])
    xticks(ticksx(ind))
    xticklabels(datestr(ticksx(ind)))
    set(gca,'TickLength',[0.005, 0.005])
    ylabel('umol/L')
    title('NSIF Nitrate','fontweight','normal','position',[datenum(year,9,7), 24, 0]);
    
    %save figure
    set(gcf,'color','w');
    export_fig('-r300',strcat(mooring_name,'_September_April_',num2str(year)),'-pdf')
end
