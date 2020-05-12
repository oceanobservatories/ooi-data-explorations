%%
% Platform:
    % Station Papa:
	% GP03FLMA, GP03FLMB, GP02HYPM
	% GPGL276, GPGL361, GPGL362, GPGL363, GPGL364, GPGL365, GPGL453, GPGL523, GPGL525, GPGL537, GPGL469
	% GPPG514, GPPG515, GPPG575, GPPG576
    %
	% Irminger Sea:
	% GI03FLMA, GI03FLMB, GI01SUMO, GI02HYPM
	% GIGL463, GIGL469, GIGL477, GIGL478, GIGL484, GIGL485, GIGL486, GIGL493, GIGL495, GIGL559, GIGL453, GIGL525, GIGL560
	% GIPG528, GIPG564, GIPG577, GIPG581
    %
    % Southern Ocean:
	% GS03FLMA, GS03FLMB, GS01SUMO, GS02HYPM
    % GSGL484, GSGL485, GSGL486, GSGL524, GSGL560, GSGL561
	% GSPG565, GSPG566
    %
    % Argentine Basin:
	% GA03FLMA, GA03FLMB, GA01SUMO, GA02HYPM
	% GAGL364, GAGL470, GAGL493, GAGL494, GAGL495, GAGL496, GAGL538
	% GAPG562, GAPG563, GAPG578, GAPG580
%Node:
    % RISER, BUOY, NSIF, PROFILER-U, PROFILER-L
%Instrument Class:
    % ADCP, CTD, DOSTA, FLORT, METBK1, METBK2, METBK1-hr, METBK2-hr, MOPAK, NUTNR, OPTAA, PARAD, PCO2A, PCO2W, PHSEN, PRESF, SPKIR, VEL3D, VELPT
	% WAVSS_Stats, WAVSS_MeanDir, WAVSS_NonDir, WAVSS_Motion, WAVSS_Fourier
%Method:
    % Telemetered, RecoveredHost, RecoveredInst, RecoveredWFP
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
%Specify metadata for Irminger Sea surface mooring
platform_name = 'GI01SUMO';
node = 'BUOY';
instrument_class = 'METBK1';
method = 'Telemetered';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);
%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);
%Get Data
[gi01sumo_metbk1_variables, mtime_metbk1, netcdfFilenames] = M2M_Data(variables, nclist);

    
%%
%Specify metadata for Irminger Sea surface mooring
platform_name = 'GI01SUMO';
node = 'BUOY';
instrument_class = 'METBK2';
method = 'Telemetered';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);
%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);
%Get Data
[gi01sumo_metbk2_variables, mtime_metbk2, netcdfFilenames] = M2M_Data(variables, nclist);

%%
%Plot data
figure('units','normalized','outerposition',[0 0 1 1])
%Plot barometric pressure
ax1 = subplot(311);
plot(mtime_metbk1,gi01sumo_metbk1_variables(7).data)
hold on
plot(mtime_metbk2,gi01sumo_metbk2_variables(7).data)
datetick('x')
title([platform_name ' ' node ' ' strrep(gi01sumo_metbk2_variables(7).name,'_',' ')])
ylabel(variables(7).units)
legend('METBK1','METBK2')

%Plot air temperature
ax2 = subplot(312);
plot(mtime_metbk1,gi01sumo_metbk1_variables(8).data)
hold on
plot(mtime_metbk2,gi01sumo_metbk2_variables(8).data)
datetick('x')
title([platform_name ' ' node ' ' strrep(gi01sumo_metbk2_variables(8).name,'_',' ')])
ylabel(variables(8).units)
legend('METBK1','METBK2')

%Plot zonal winds
ax3 = subplot(313);
plot(mtime_metbk1,gi01sumo_metbk1_variables(5).data)
hold on
plot(mtime_metbk2,gi01sumo_metbk2_variables(5).data)
datetick('x')
title([platform_name ' ' node ' ' strrep(gi01sumo_metbk2_variables(5).name,'_',' ')])
ylabel(variables(5).units)
legend('METBK1','METBK2')
linkaxes([ax1 ax2 ax3],'x')
