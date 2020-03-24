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
start_date='2019-10-01T00:00:00.000Z';
end_date='2019-11-10T23:59:59.000Z';

%%
%Specify metadata
platform_name = 'CE09OSSM';
node = 'MFN';
instrument_class = 'OPTAA';
method = 'Telemetered';

%Get M2M URL
[uframe_dataset_name,variables]=M2M_URLs(platform_name,node,instrument_class,method);
%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);
%Get Data
[~, ~, netcdfFilenames] = M2M_Data(variables, nclist, false);

%Process OPTAA Data
tf_addDiscreteWavelengthTimeSeries = true;
optaa_variables = M2M_OPTAA(netcdfFilenames, tf_addDiscreteWavelengthTimeSeries);

%Average all spectra data
attenuation_spectra_1d_mean = nanmean(optaa_variables.attenuation_spectra_2d,2);
absorption_spectra_1d_mean = nanmean(optaa_variables.absorption_spectra_2d,2);

%Calculate spectra stddev
attenuation_spectra_1d_std = nanstd(optaa_variables.attenuation_spectra_2d,0,2);
absorption_spectra_1d_std = nanstd(optaa_variables.absorption_spectra_2d,0,2);

%Plot the data
mtime      = optaa_variables.mtime;
wavelength = optaa_variables.wavelength;
figure('units','normalized','outerposition',[0 0 1 1])
subplot(221)
pcolor(mtime,wavelength,optaa_variables.attenuation_spectra_2d)
shading flat
caxis([0 .4])
colormap(jet)
c=colorbar;
title(c,'m^-^1')
ylabel('Wavelength (nm)')
title({[platform_name ' ' node ' ' 'Attenuation'],'pre-deployment calibration offsets have not been applied to these data'})
datetick('x',1,'keeplimits')
set(gca,'Layer','top')
%
subplot(222)
ind=~isnan(attenuation_spectra_1d_mean);
upper = attenuation_spectra_1d_mean(ind)'+attenuation_spectra_1d_std(ind);
lower = attenuation_spectra_1d_mean(ind)'-attenuation_spectra_1d_std(ind);
fill([wavelength(ind)' fliplr(wavelength(ind)')], [upper fliplr(lower)],...
    [.9 .9 .9], 'linestyle', 'none')
hold on
plot(wavelength(ind),attenuation_spectra_1d_mean(ind),'k')
xlabel('Wavelength (nm)')
ylabel('Attenuation (m^-^1)')
xlim([400 744])
ylim([0 .4])
title(['Time Period: ' datestr(min(mtime)) ' to ' datestr(max(mtime))])
%
subplot(223)
pcolor(mtime,wavelength,optaa_variables.absorption_spectra_2d)
shading flat
caxis([0 .4])
colormap(jet)
c=colorbar;
title(c,'m^-^1')
ylabel('Wavelength (nm)')
title({[platform_name ' ' node ' ' 'Absorption'],'pre-deployment calibration offsets have not been applied to these data'})
datetick('x',1,'keeplimits')
set(gca,'Layer','top')
%
subplot(224)
ind=~isnan(absorption_spectra_1d_mean);
upper = absorption_spectra_1d_mean(ind)'+absorption_spectra_1d_std(ind);
lower = absorption_spectra_1d_mean(ind)'-absorption_spectra_1d_std(ind);
fill([wavelength(ind)' fliplr(wavelength(ind)')], [upper fliplr(lower)],...
    [.9 .9 .9], 'linestyle', 'none')
hold on
plot(wavelength(ind),absorption_spectra_1d_mean(ind),'k')
xlabel('Wavelength (nm)')
ylabel('Absorption (m^-^1)')
xlim([400 744])
ylim([0 .4])
title(['Time Period: ' datestr(min(mtime)) ' to ' datestr(max(mtime))])

%.. plot time series of individual channels
%.. .. the 1D channel wavelengths are interpolated from the 2D data and are
%.. .. hardcoded in the "revamp" function (and can be changed by the user)
if tf_addDiscreteWavelengthTimeSeries
    figure('units','normalized','outerposition',[0 0 1 1])
    %.. get the 1D channel data by converting the structure to a cell, then
    %.. converting data in the cell selected using fields to numeric
    celldata = struct2cell(optaa_variables);
    fields = fieldnames(optaa_variables);
   %.. absorption
    tf_abs = contains(fields, 'abs_');
    absChannels = cell2mat(celldata(tf_abs));
    %.. beam attenuation
    tf_ccc = contains(fields, 'beam_c_');
    cccChannels = cell2mat(celldata(tf_ccc));
    %.. get 1D wavelengths
    waveChar = char(fields(tf_abs));
    waveChar(:, 8) = ' ';
    waveChar = waveChar(:, 5:8);
    waveNumeric = str2num(waveChar);  %#ok   column vector of numbers
    waveChar = waveChar';
    waveChar = waveChar(:)';  % row character vector
    
    subplot(211)
    plot(mtime, absChannels)
    title({[platform_name ' ' node ' ' 'Absorption Wavelengths'] ... 
        'pre-deployment calibration offsets have not been applied to these data'})
    ylabel('Absorption (m^-^1)')
    legend(split(waveChar(1:end-1),' '))
    datetick('x',1,'keeplimits')
    subplot(212)
    plot(mtime, cccChannels)
    datetick('x',1,'keeplimits')
    ylabel('Absorption (m^-^1)')
    legend(split(waveChar(1:end-1),' '))
    title({[platform_name ' ' node ' ' 'Beam Attenuation Wavelengths'] ...
        'pre-deployment calibration offsets have not been applied to these data'})
end
