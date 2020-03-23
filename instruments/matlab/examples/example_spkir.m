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
close all
clearvars

%.. set login and URL details
username = 'OOIAPI-853A3LA6QI3L62';  %api_key
password = 'WYAN89W5X4Z0QZ';	%api_token
% Set Authorization header field in weboptions
options = weboptions('CertificateFilename','','HeaderFields',{'Authorization',...
    ['Basic ' matlab.net.base64encode([username ':' password])]}, 'Timeout', 120);

%.. set time period of interest
start_date='2019-03-01T00:00:00.000Z';
end_date='2019-05-31T23:59:59.000Z';

%%
%Specify metadata
mooring_name = 'CE02SHSM';
node = 'NSIF';
instrument_class = 'SPKIR';
method = 'RecoveredHost';

%Get M2M URL
[uframe_dataset_name,variables]=M2M_URLs(mooring_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
[~, ~, netcdfFilenames] = M2M_Data(variables, nclist, false);

%%
%Process SPKIR Data
[spkir_variables] = M2M_SPKIR(netcdfFilenames);

%Plot the data
%.. time series of all wavelength channels
figure('units','normalized','outerposition',[0 0 1 1])
subplot(211)
plot(spkir_variables.mtime_binned,spkir_variables.Ed_412nm_binned,'color',[85/255 26/255 139/255])
hold on
plot(spkir_variables.mtime_binned,spkir_variables.Ed_443nm_binned,'color',[0/255 0/255 255/255])
plot(spkir_variables.mtime_binned,spkir_variables.Ed_490nm_binned,'color',[0/255 139/255 139/255])
plot(spkir_variables.mtime_binned,spkir_variables.Ed_510nm_binned,'color',[80/255 200/255 120/255])
plot(spkir_variables.mtime_binned,spkir_variables.Ed_555nm_binned,'color',[173/255 255/255 47/255])
plot(spkir_variables.mtime_binned,spkir_variables.Ed_620nm_binned,'color',[255/255 165/255 0/255])
plot(spkir_variables.mtime_binned,spkir_variables.Ed_683nm_binned,'color',[255/255 0/255 0/255])
legend('412nm','443nm','490nm','510nm','555nm','620nm','683nm')
ylim([0 100])
datetick('x',1,'keeplimits')
ylabel('uW cm-2 nm-1')
title([mooring_name ' ' node ' ' 'Downwelling Spectral Irradiance;  Local time is UTC-7'])
%.. spectrum at time of maximum measured irradiance

subplot(212)
%Plot spectrum when highest downwelling irradiance observed; 
celldata = struct2cell(spkir_variables);
%.. last 7 fields contain the binned data
celldata = celldata(end-6:end);
%.. time series data must be contained in row vectors:
binnedChannelArray = cell2mat(celldata);
[~, idxMax] = max(sum(binnedChannelArray, 1));
timeAtMax     = spkir_variables.mtime_binned(idxMax);
spectrumMax = binnedChannelArray(:, idxMax);  % column vector
plot(spkir_variables.wavelength, spectrumMax, 'k*-')
xlabel('Wavelength (nm)')
ylabel('uW cm-2 nm-1')
title({[mooring_name ' ' node],['Maximum brightness at 7m depth was observed on ' datestr(timeAtMax) ' UTC;  Local time is UTC-7']})
