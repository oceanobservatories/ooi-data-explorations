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
    % GSGL486, GSGL560
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
start_date='2015-09-01T00:00:00.000Z';
end_date='2015-10-31T23:59:59.000Z';

%%
%Specify metadata
platform_name = 'GSGL486';
node = 'GLIDER';
instrument_class = 'CTD';
method = 'RecoveredHost';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);

%Make M2M Call
[nclist] = M2M_Call(uframe_dataset_name,start_date,end_date,options);

%Get Data
%[ctd_variables, ctd_mtime, netcdfFilenames] = M2M_Data(variables, nclist, false);   %This will download .nc file(s) and read in the data from the local files
[ctd_variables, ctd_mtime, netcdfFilenames] = M2M_Data(variables, nclist);  %This will use the opendap to read in the data from remote files

%Example plot
figure
scatter(ctd_mtime,ctd_variables(5).data,5,ctd_variables(2).data,'filled')
caxis([4 6])
c=colorbar;
title(c,ctd_variables(2).units)
set(gca, 'YDir','reverse')
ylabel(ctd_variables(5).units)
ylim([-10 1000])
datetick('x',1)
title([platform_name ' ' strrep(ctd_variables(2).name,'_',' ')])
box on

figure
scatter(ctd_variables(8).data,ctd_variables(7).data,4,ctd_mtime,'filled')
ylabel([ctd_variables(7).name])
xlabel([ctd_variables(8).name])
c=colorbar;
title(c,'time')
title([platform_name ' ' start_date ' ' '--' ' ' end_date])
box on
