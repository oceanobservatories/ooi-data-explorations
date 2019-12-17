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
start_date='2018-12-01T00:00:00.000Z';
end_date='2018-12-31T23:59:59.000Z';

%%
%Specify metadata
mooring_name = 'CE02SHSM';
node = 'BUOY';
instrument_class = 'WAVSS_MeanDir';
method = 'RecoveredHost';

%Get M2M URL
[uframe_dataset_name,variables]=M2M_URLs(mooring_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
%[wavss_variables, wavss_mtime, netcdfFilenames] = M2M_Data(variables, nclist, false);   %This will download .nc file(s) and read in the data from the local files
[wavss_variables, wavss_mtime, netcdfFilenames] = M2M_Data(variables, nclist);  %This will use the opendap to read in the data from remote files

[row,col]=size(wavss_variables(12).data);  %wavss_a_corrected_directional_wave_direction
frequency_matrix = nan(row,col);
initial_frequency = wavss_variables(4).data;  %initial_frequency
frequency_interval = wavss_variables(5).data; %frequency_spacing

for i = 1:col
    k=-1;
    for j = 1:row
        k=k+1;
        frequency_matrix(j,i)=initial_frequency(i)+(frequency_interval(i)*k);
    end    
end

%Example plots
figure
pcolor(wavss_mtime,frequency_matrix(:,1),wavss_variables(12).data); %wavss_a_corrected_directional_wave_direction
shading flat
ylim([0.05 .6])
ylabel('Hz')
datetick('x')
title(strrep(wavss_variables(12).name,'_','-'))
c=colorbar;
title(c,wavss_variables(12).units)

figure
pcolor(wavss_mtime,frequency_matrix(:,1),wavss_variables(6).data); %psd_mean_directional
shading flat
ylim([0.05 .6])
ylabel('Hz')
datetick('x')
caxis([0 1.2])
title(strrep(wavss_variables(6).name,'_','-'))
c=colorbar;
title(c,wavss_variables(6).units)

%now bin energy by direction and frequency to produce polar plot 
dir_bins = 0:10:370;
hz_bins = frequency_matrix(:,1)';
len= size(frequency_matrix,1)*size(frequency_matrix,2);
dir=reshape(wavss_variables(12).data,[len 1]); %wavss_a_corrected_directional_wave_direction
e=  reshape(wavss_variables(6).data, [len 1]); %psd_mean_directional
hz= reshape(frequency_matrix,        [len 1]);

energy_vs_direction=zeros(length(hz_bins)-1,length(dir_bins)-1);
for ii = 1:length(dir_bins)-1
    for jj = 1:length(hz_bins)-1
        ind = find(dir>=dir_bins(ii) & dir<dir_bins(ii+1) & hz>=hz_bins(jj) & hz<hz_bins(jj+1));
        if ~isempty(ind) 
            energy_vs_direction(jj,ii)=nanmean(e(ind));
        end
    end
end

figure;
pcolor(dir_bins(1:37),hz_bins(1:122),energy_vs_direction);
caxis([0 1]);
shading flat
colorbar; ylabel('Hz');
set(gca,'xtick',0:90:360); xlabel('heading');  

figure; title('wave rose');
polarPcolor(hz_bins(1:122)+0.0075,dir_bins(1:37),energy_vs_direction);

PolarContour(energy_vs_direction(1:60,:)', hz_bins(1:60)')
caxis([0 12])
c=colorbar;title(c,'m^2/Hz')
