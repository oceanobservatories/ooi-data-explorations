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
instrument_class = 'FLORT';
method = 'RecoveredHost';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);

for ii = 1:2
    %Make M2M Call
    [nclist] = M2M_Call(uframe_dataset_name(ii),start_date,end_date,options);
    %Get Data
    [flort_variables{ii}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);  %change nclist(1) to nclist to get it to fail
end

%Example plot
figure
subplot(211)
scatter((flort_variables{1,1}(1).data)/60/60/24+datenum(1900,1,1),...
        (flort_variables{1,1}(13).data),4,flort_variables{1,1}(3).data,'filled')
caxis([0 1])
c=colorbar;
title(c,flort_variables{1,1}(3).units)
set(gca, 'YDir','reverse')
ylabel(flort_variables{1,1}(13).units)
ylim([-2 210])
datetick('x',1)
title([platform_name ' ' strrep(flort_variables{1,1}(3).name,'_',' ')])
box on

subplot(212)
scatter((flort_variables{1,1}(1).data)/60/60/24+datenum(1900,1,1),...
        (flort_variables{1,1}(13).data),4,flort_variables{1,1}(6).data,'filled')
caxis([0 6])
c=colorbar;
title(c,flort_variables{1,1}(6).units)
set(gca, 'YDir','reverse')
ylabel(flort_variables{1,1}(13).units)
ylim([-2 210])
datetick('x',1)
title([platform_name ' ' strrep(flort_variables{1,1}(6).name,'_',' ')])
box on

figure
scatter(flort_variables{1,1}(15).data,flort_variables{1,1}(14).data,4,...
    (flort_variables{1,1}(1).data)/60/60/24+datenum(1900,1,1),'filled')
ylabel([flort_variables{1,1}(14).name])
xlabel([flort_variables{1,1}(15).name])
c=colorbar;
title(c,'time')
title([platform_name ' ' start_date ' ' '--' ' ' end_date])
box on
