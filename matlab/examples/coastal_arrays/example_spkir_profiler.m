%.. 2019-12-18: CMRisien. original code: example_spkir.m
%.. 2020-06-15: RADesiderio. added logic so that calls to profiler data are not median filtered.
%.. 2020-06-19: RADesiderio. changed plotting and renamed to example_spkir_profiler.m
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
start_date = '2019-07-01T00:00:00.000Z';
end_date   = '2019-07-14T23:59:59.999Z';

%%
%Specify metadata
platform_name    = 'CE01ISSP';
node             = 'PROFILER';
instrument_class = 'SPKIR';
method           = 'RecoveredCSPP';

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
%Plot the profiler data
vars = fieldnames(spkir_variables);
if ~any(contains(vars, 'pressure'))
    tf_noPressureData = true;
else
    tf_noPressureData = false;
end

siteDepth = 25;  % for removing spike artifacts in pressure

color{1} = [144/255    0/255  255/255];
color{2} = [  0/255    0/255  255/255];
color{3} = [  0/255  155/255  175/255];
color{4} = [ 80/255  200/255  120/255];
color{5} = [ 96/255  244/255   48/255];
color{6} = [255/255  165/255    0/255];
color{7} = [255/255    0/255    0/255];

%.. transfer irradiance data from structure fields into matrix
tf_Ed            = contains(vars, 'Ed_');
celldata         = struct2cell(spkir_variables);
celldata(~tf_Ed) = [];
irradiance       = cell2mat(celldata)';
%.. 7 columns of data, one for each wavelength

%.. 1st plot: unfiltered time series of pressure and all wavelength channels
figure('units','normalized','outerposition',[0 0 1 1])
if tf_noPressureData
    subplot(421);
    ylabel('pressure [db]');
    title('NO PRESSURE DATA');
else
    ax(8) = subplot(421);  % upper left subplot of 4x2 grid
    plot(spkir_variables.mtime, spkir_variables.pressure, 'kx-')
    axis ij
    datetick('x',1,'keeplimits')
    title([instrument_class ',  ' platform_name ' ' node ': Local time is UTC-7'])
    ylabel('pressure [db]');
end

waveString = string(num2str(spkir_variables.wavelength));
%.. order subplots down the 1st column then down the second
orderSubplots = [3 5 7 2 4 6 8]; 
for ii = 1:7
    ax(ii) = subplot(4, 2, orderSubplots(ii) );
    plot(spkir_variables.mtime, irradiance(:, ii), 'color', color{ii}, 'Marker', 'x')
    datetick('x',1,'keeplimits')
    ylabel('uW cm-2 nm-1')
    title(['Downwelling Irradiance at ' char(waveString(ii)) 'nm'], 'Color', color{ii});
end
linkaxes(ax, 'x');

%.. 2nd plot: plots of profiles for all wavelength channels
if tf_noPressureData
    disp('2nd plot cannot be made because there are no pressure data.');
    return
end
pressure = spkir_variables.pressure;
pressure(pressure > siteDepth + 5) = nan;
%.. set up irradiance data and limits for log plotting
%.. .. measurements in the absence of light can result in small negative values.
%.. .. values below 1 may be significant in certain cases but are disregarded for log plots.
irradianceMinLim = 1;  % uW/cm^2/nm
irradiance(irradiance <= irradianceMinLim) = 1;
irradiance = log10(irradiance);
irr_logMax = nanmax(irradiance, [], 1);

figure('units','normalized','outerposition',[0 0 1 1])
axx(8) = subplot(421);
plot(spkir_variables.mtime, pressure, 'kx-')
axis ij
datetick('x',1,'keeplimits')
title([instrument_class ',  ' platform_name ' ' node ': Local time is UTC-7'])
ylabel('pressure [db]');
for ii = 1:7
    axx(ii) = subplot(4, 2, orderSubplots(ii) );
    fastscatter(spkir_variables.mtime', pressure, irradiance(:, ii), 'MarkerSize', 16)
    axis ij
    datetick('x',1,'keeplimits')
    ylabel('pressure [db]')
    colorbar
    caxis([0 irr_logMax(ii)]);
    title(['log10(Downwelling Irradiance) [uW cm-2 nm-1] at ' char(waveString(ii)) 'nm'], 'Color', color{ii});
end
linkaxes(axx, 'x');
