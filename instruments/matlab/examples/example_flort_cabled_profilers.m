%%
% Mooring/Glider:
    %CE01ISSM, CE01ISSP, CE02SHSM, CE02SHBP, CE02SHSP, CE04OSSM, CE04OSBP, CE06ISSM, CE06ISSP, CE07SHSM, CE07SHSP, CE09OSSM, CE09OSPM
    %CEGL386, CEGL384, CEGL383, CEGL382, CEGL381, CEGL327, CEGL326, CEGL320, CEGL319, CEGL312, CEGL311, CEGL247
%Node:
    %BUOY, NSIF, MFN, BEP, PROFILER, GLIDER
%Instrument Class:
    %ADCP, CTD, DOSTA, FDCHP, FLORT, METBK1, METBK2, METBK1-hr, METBK2-hr, MOPAK, NUTNR, OPTAA, PARAD, PCO2A, PCO2W, PHSEN, PRESF, SPKIR, VEL3D, VELPT, ZPLSC
	%WAVSS_Stats, WAVSS_MeanDir, WAVSS_NonDir, WAVSS_Motion, WAVSS_Fourier
%Method:
    %Telemetered, RecoveredHost, RecoveredInst, RecoveredCSPP, RecoveredWFP, Streamed
%%
close all
clearvars

%.. set login and URL details
username = 'OOIAPI-853A3LA6QI3L62';  %api_key
password = 'WYAN89W5X4Z0QZ';	%api_token
% Set Authorization header field in weboptions
options = weboptions('CertificateFilename','','HeaderFields',{'Authorization',...
    ['Basic ' matlab.net.base64encode([username ':' password])]}, 'Timeout', 120);

%.. set time period of interest
start_date='2019-07-09T00:00:00.000Z';
end_date='2019-07-10T23:59:59.000Z';

%%
%Specify metadata
mooring_name = 'CE04OSPD';
node = 'PROFILER';
instrument_class = 'FLORT';
method = 'RecoveredInst';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(mooring_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
[CE04OSPD_variables, CE04OSPD_mtime, netcdfFilenames] = M2M_Data(variables, nclist);  %This will use the opendap to read in the data from remote files

%Example plot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(211)
plot(CE04OSPD_mtime,CE04OSPD_variables(2).data)
datetick('x',1)
title([mooring_name ' ' node ' ' strrep(CE04OSPD_variables(2).name,'_',' ')])
ylabel(CE04OSPD_variables(2).units)

%%
%Specify metadata
mooring_name = 'CE04OSPS';
node = 'PROFILER';
instrument_class = 'FLORT';
method = 'Streamed';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(mooring_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
[CE04OSPS_variables, CE04OSPS_mtime, netcdfFilenames] = M2M_Data(variables, nclist);  %This will use the opendap to read in the data from remote files

%Example plot
figure(1)
subplot(212)
scatter(CE04OSPS_mtime,CE04OSPS_variables(7).data,8,CE04OSPS_variables(3).data)
caxis([0 1])
cb=colorbar;title(cb,CE04OSPS_variables(3).units)
set(cb,'Position',[.92 .18 0.01 0.2])
set(gca, 'YDir','reverse')
datetick('x',1)
title([mooring_name ' ' node ' ' strrep(CE04OSPS_variables(3).name,'_',' ')])
ylabel(CE04OSPS_variables(7).units)
box on
