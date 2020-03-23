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
start_date='2016-06-01T00:00:00.000Z';
end_date='2016-06-21T23:59:59.000Z';

%%
%Specify metadata
mooring_name = 'CE07SHSM';
node = 'BUOY';
instrument_class = 'METBK1';
method = 'RecoveredHost';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(mooring_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
%[metbk_variables, metbk_mtime, netcdfFilenames] = M2M_Data(variables, nclist, false);   %This will download .nc file(s) and read in the data from the local files
[metbk_variables, metbk_mtime, netcdfFilenames] = M2M_Data(variables, nclist);  %This will use the opendap to read in the data from remote files

%%
%Example plot
figure(1)
%Plot W/E Winds
subplot(311)
plot(metbk_mtime,metbk_variables(5).data)
xlim([min(metbk_mtime) max(metbk_mtime)])
datetick('x',1,'keeplimits')
title([mooring_name ' ' node ' ' strrep(metbk_variables(5).name,'_',' ')])
ylabel(metbk_variables(5).units)
ylim([-10 10])
%Plot N/S Winds
subplot(312)
plot(metbk_mtime,metbk_variables(6).data)
xlim([min(metbk_mtime) max(metbk_mtime)])
datetick('x',1,'keeplimits')
title([mooring_name ' ' node ' ' strrep(metbk_variables(6).name,'_',' ')])
ylabel(metbk_variables(6).units)
ylim([-10 10])
%Plot Barometric Pressure
subplot(313)
plot(metbk_mtime,metbk_variables(7).data)
xlim([min(metbk_mtime) max(metbk_mtime)])
datetick('x',1,'keeplimits')
title([mooring_name ' ' node ' ' strrep(metbk_variables(7).name,'_',' ')])
ylabel(metbk_variables(7).units)

%Plot a wind stickplot
figure(2)
subplot(211)
h=stickplot(metbk_mtime(1:10:end),metbk_variables(5).data(1:10:end),metbk_variables(6).data(1:10:end),[min(metbk_mtime) max(metbk_mtime) -10 10]);
set(h,'color','k','linewidth',.4)
ylabel(metbk_variables(6).units)
datetick('x',1,'keeplimits')
box on
title([mooring_name ' ' node ' Winds'])
