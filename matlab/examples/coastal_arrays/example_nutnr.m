%%
% Platform:
    % Endurance Array:
	% CE01ISSM, CE01ISSP, CE02SHSM, CE02SHBP, CE02SHSP, CE04OSSM, CE04OSBP, CE06ISSM, CE06ISSP, CE07SHSM, CE07SHSP, CE09OSSM, CE09OSPM
    % CEGL386, CEGL384, CEGL383, CEGL382, CEGL381, CEGL327, CEGL326, CEGL320, CEGL319, CEGL312, CEGL311, CEGL247
	%
	% Pioneer Array:
	% CP01CNPM, CP01CNSM, CP02PMCI, CP02PMCO, CP02PMUI, CP02PMUO, CP03ISPM, CP03ISSM, CP04OSPM, CP04OSSM
	% CPGL335, CPGL336, CPGL339, CPGL340, CPGL374, CPGL375, CPGL376, CPGL379, CPGL380, CPGL387, CPGL388, CPGL389, CPGL514
%Node:
    % BUOY, NSIF, MFN, BEP, PROFILER, GLIDER
%Instrument Class:
    % ADCP, CTD, DOSTA, FDCHP, FLORT, METBK1, METBK2, METBK1-hr, METBK2-hr, MOPAK, NUTNR, OPTAA, PARAD, PCO2A, PCO2W, PHSEN, PRESF, SPKIR, VEL3D, VELPT, ZPLSC
	% WAVSS_Stats, WAVSS_MeanDir, WAVSS_NonDir, WAVSS_Motion, WAVSS_Fourier
%Method:
    % Telemetered, RecoveredHost, RecoveredInst, RecoveredCSPP, RecoveredWFP, Streamed
%%
% load the access credentials
try
    load('ooinet.credentials.mat')  % returns a variable called options
catch
    error(['Unable to load access credentials. Users need to create a ' ...
           'weboptions object with their personal OOINet API keys. See ' ...
           'README for more information on how to create this.'])
end %try

%.. set time period of interest
start_date='2015-04-01T00:00:00.000Z';
end_date='2019-09-30T23:59:59.000Z';
% start_date = '2016-05-16T00:00:00.000Z';
% end_date   = '2019-06-22T23:59:59.000Z';

%%
%Specify metadata
platform_name = 'CE02SHSM';
node = 'NSIF';
instrument_class = 'NUTNR';
method = 'Telemetered';
%method = 'RecoveredHost';
%method = 'RecoveredInst';

%Get M2M URL
[uframe_dataset_name,variables]=M2M_URLs(platform_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
%[variables, mtime, netcdfFilenames] = M2M_Data(variables, nclist, false);   %This will download .nc file(s) and read in the data from the local files
[variables, mtime, netcdfFilenames] = M2M_Data(variables, nclist);  %This will use the opendap to read in the data from remote files

%%
%Example plot
if strcmpi(method, 'Telemetered')
    figure(101)
elseif strcmpi(method, 'RecoveredHost')
    figure(102)
elseif strcmpi(method, 'RecoveredInst')
    figure(103)
end

subplot(211)
plot(mtime,variables(2).data, 'x-')
datetick('x',1)
title(['L1: ' strrep(variables(2).name,'_',' ')])
subplot(212)
plot(mtime,variables(3).data, 'x-')
datetick('x',1)
title(['L2: ' strrep(variables(3).name,'_',' ')])
linkaxes([subplot(211) subplot(212)], 'x')
