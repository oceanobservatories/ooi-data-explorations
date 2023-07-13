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
% load the access credentials
try
    load('ooinet.credentials.mat')  % returns a variable called options
catch
    error(['Unable to load access credentials. Users need to create a ' ...
           'weboptions object with their personal OOINet API keys. See ' ...
           'README for more information on how to create this.'])
end %try

%.. set time period of interest
start_date='2016-06-01T00:00:00.000Z';
end_date='2016-07-31T23:59:59.000Z';

%%
%Specify metadata
platform_name = 'GSPG566';
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
subplot(211)
scatter(ctd_mtime,ctd_variables(5).data,5,ctd_variables(2).data,'filled')
caxis([5 7])
c=colorbar;
title(c,ctd_variables(2).units)
set(gca, 'YDir','reverse')
ylabel(ctd_variables(5).units)
ylim([-2 210])
datetick('x',1)
title([platform_name ' ' strrep(ctd_variables(2).name,'_',' ')])
box on
subplot(212)
scatter(ctd_mtime,ctd_variables(5).data,5,ctd_variables(3).data,'filled')
caxis([34 34.4])
c=colorbar;
title(c,ctd_variables(3).units)
set(gca, 'YDir','reverse')
ylabel(ctd_variables(5).units)
ylim([-2 210])
datetick('x',1)
title([platform_name ' ' strrep(ctd_variables(3).name,'_',' ')])
box on

figure
scatter(ctd_variables(8).data,ctd_variables(7).data,4,ctd_mtime,'filled')
ylabel([ctd_variables(7).name])
xlabel([ctd_variables(8).name])
c=colorbar;
title(c,'time')
title([platform_name ' ' start_date ' ' '--' ' ' end_date])
box on
