%.. 2019-12-18: CMRisien. original code: example_spkir.m
%.. 2020-06-15: RADesiderio. added logic so that calls to profiler data are not median filtered.
%.. this function will be renamed to example_spkir_mooring in a future pull request
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
start_date = '2019-03-01T00:00:00.000Z';
end_date   = '2019-05-31T23:59:59.999Z';

%%
%Specify metadata
platform_name    = 'CE02SHSM';
node             = 'NSIF';
instrument_class = 'SPKIR';
method           = 'RecoveredHost';

%Get M2M URL
[uframe_dataset_name,variables]=M2M_URLs(platform_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
%.. download netcdf files (3rd argument false) for M2M_SPKIR
[~, ~, netcdfFilenames] = M2M_Data(variables, nclist, false);

%Process SPKIR Data
if strcmpi(node, 'PROFILER')
    tf_medianFilter = false;
else
    tf_medianFilter = true;
end
[spkir_variables] = M2M_SPKIR(netcdfFilenames, tf_medianFilter);

%%
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
title([platform_name ' ' node ' ' 'Downwelling Spectral Irradiance;  Local time is UTC-7'])
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
title({[platform_name ' ' node],['Maximum brightness at 7m depth was observed on ' datestr(timeAtMax) ' UTC;  Local time is UTC-7']})
