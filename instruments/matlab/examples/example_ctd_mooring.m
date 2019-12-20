%%
% Mooring/Glider:
    %CE01ISSM, CE01ISSP, CE02SHSM, CE02SHBP, CE02SHSP, CE04OSSM, CE04OSBP, CE06ISSM, CE06ISSP, CE07SHSM, CE07SHSP, CE09OSSM, CE09OSPM
    %CEGL386, CEGL384, CEGL383, CEGL382, CEGL381, CEGL327, CEGL326, CEGL320, CEGL319, CEGL312, CEGL311, CEGL247
%Node:
    %BUOY, NSIF, MFN, BEP, PROFILER, GLIDER
%Instrument Class:
    %ADCP, CTD, DOSTA, FDCHP, FLORT, METBK, METBK-hr, MOPAK, NUTNR, OPTAA, PARAD, PCO2A, PCO2W, PHSEN, PRESF, SPKIR, VEL3D, VELPT, ZPLSC
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
start_date='2019-01-01T00:00:00.000Z';
end_date='2019-09-30T23:59:59.000Z';

%%
%Specify metadata
mooring_name = 'CE07SHSM';
node = 'NSIF';
instrument_class = 'CTD';
method = 'RecoveredInst';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(mooring_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
%[variables, mtime, netcdfFilenames] = M2M_Data(variables, nclist, false);   %This will download .nc file(s) and read in the data from the local files
[variables, mtime, netcdfFilenames] = M2M_Data(variables, nclist);  %This will use the opendap to read in the data from remote files

%Example plot
plot(mtime,variables(2).data)
datetick('x',1)
title([mooring_name ' ' node ' ' strrep(variables(2).name,'_',' ')])
ylabel(variables(2).units)
