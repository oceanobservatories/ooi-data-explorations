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
start_date='2015-01-01T00:00:00.000Z';
end_date='2016-12-30T23:59:59.000Z';

%%
%Specify metadata for Irminger Sea Flanking Mooring A
platform_name = 'GI03FLMA';
node = 'RISER';
instrument_class = 'VELPT';
method = 'RecoveredInst';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);

for ii = 1:2
    %Make M2M Call
    [nclist] = M2M_Call(uframe_dataset_name{ii},start_date,end_date,options);
    %Get Data
    [gi03flma_velpt_variables{ii}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
end

%%
%Specify metadata for Irminger Sea Flanking Mooring B
platform_name = 'GI03FLMB';
node = 'RISER';
instrument_class = 'VELPT';
method = 'RecoveredInst';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);

for ii = 1:4
    %Make M2M Call
    [nclist] = M2M_Call(uframe_dataset_name{ii},start_date,end_date,options);
    %Get Data
    [gi03flmb_velpt_variables{ii}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
end

%%
%Plot data
figure('units','normalized','outerposition',[0 0 1 1])
subplot(211)
hold on
for i = 1:4
    scatter((gi03flma_velpt_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        (gi03flma_velpt_variables{1,i}(9).data)*0.001,4,gi03flma_velpt_variables{1,i}(3).data,'filled')
end
datetick('x',1,'keeplimits')
ylabel('dbar')
ylim([1600 2800])
set(gca, 'YDir','reverse')
caxis([-.4 .4])
c=colorbar;
title(c,gi03flma_velpt_variables{1,i}(3).units)
box on
title('GI03FLMA Meridional Velocities','fontweight','normal')

subplot(212)
hold on
for i = 1:4
    scatter((gi03flmb_velpt_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        (gi03flmb_velpt_variables{1,i}(9).data)*0.001,4,gi03flmb_velpt_variables{1,i}(3).data,'filled')
end
datetick('x',1,'keeplimits')
ylabel('dbar')
ylim([1600 2800])
set(gca, 'YDir','reverse')
caxis([-.4 .4])
c=colorbar;
title(c,gi03flmb_velpt_variables{1,i}(3).units)
box on
title('GI03FLMB Meridional Velocities','fontweight','normal')
