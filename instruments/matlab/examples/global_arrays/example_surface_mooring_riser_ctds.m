%%
% Platform:
    % Station Papa:
	% GP03FLMA, GP03FLMB, GP02HYPM
    %
	% Irminger Sea:
	% GI03FLMA, GI03FLMB, GI01SUMO, GI02HYPM
    %
    % Southern Ocean:
	% GS03FLMA, GS03FLMB, GS01SUMO, GS02HYPM
    %
    % Argentine Basin:
	% GA03FLMA, GA03FLMB, GA01SUMO, GA02HYPM
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
start_date='2016-01-01T00:00:00.000Z';
end_date='2018-12-30T23:59:59.000Z';

%%
%Specify metadata for Irminger Sea surface mooring 
platform_name = 'GI01SUMO';
node = 'RISER';
instrument_class = 'CTD';
method = 'RecoveredInst';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);

for ii = 1:13
    disp(ii)
    %Make M2M Call
    [nclist] = M2M_Call(uframe_dataset_name{ii},start_date,end_date,options);
    %Get Data
    [gi01sumo_ctd_variables{ii}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
end

%%
%Plot seawater temperature data
figure('units','normalized','outerposition',[0 0 1 1])
subplot(211)
hold on
for i = 1:13
    scatter((gi01sumo_ctd_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        gi01sumo_ctd_variables{1,i}(6).data,4,gi01sumo_ctd_variables{1,i}(2).data,'filled')
    hold on
    scatter((gi01sumo_ctd_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        gi01sumo_ctd_variables{1,i}(7).data,4,gi01sumo_ctd_variables{1,i}(3).data,'filled')
end
datetick('x',1,'keeplimits')
ylabel(gi01sumo_ctd_variables{1,i}(6).units)
ylim([0 3000])
set(gca, 'YDir','reverse')
caxis([0 12])
c=colorbar;
title(c,gi01sumo_ctd_variables{1,i}(2).units)
box on
title('GI01SUMO Seawater Temperature','fontweight','normal')

subplot(212)
hold on
for i = 1:13
    scatter((gi01sumo_ctd_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        gi01sumo_ctd_variables{1,i}(6).data,4,gi01sumo_ctd_variables{1,i}(4).data,'filled')
    hold on
    scatter((gi01sumo_ctd_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        gi01sumo_ctd_variables{1,i}(7).data,4,gi01sumo_ctd_variables{1,i}(4).data,'filled')
end
datetick('x',1,'keeplimits')
ylabel(gi01sumo_ctd_variables{1,i}(5).units)
ylim([0 3000])
set(gca, 'YDir','reverse')
caxis([34.8 35.2])
c=colorbar;
title(c,gi01sumo_ctd_variables{1,i}(4).units)
box on
title('GI01SUMO Seawater Salinity','fontweight','normal')
