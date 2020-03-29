%%
% Platform:
    % Station Papa:
	% GP03FLMA, GP03FLMB
    %
	% Irminger Sea:
	% GI03FLMA, GI03FLMB
    %
    % Southern Ocean:
	% GS03FLMA, GS03FLMB
    %
    % Argentine Basin:
	% GA03FLMA, GA03FLMB
%Node:
    % RISER
%Instrument Class:
    % ADCP, CTD, DOSTA, FLORT, PHSEN, VELPT
%Method:
    % Telemetered, RecoveredHost, RecoveredInst
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

for ii = 1:4
    %Make M2M Call
    [nclist] = M2M_Call(uframe_dataset_name{ii},start_date,end_date,options);
    %Get Data
    [gs03flma_velpt_variables{ii}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
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
    [gs03flmb_velpt_variables{ii}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(211)
hold on
for i = 1:4
    scatter((gs03flma_velpt_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        (gs03flma_velpt_variables{1,i}(9).data)*0.001,4,gs03flma_velpt_variables{1,i}(3).data,'filled')
end
datetick('x',1,'keeplimits')
ylabel('dbar')
ylim([1600 2800])
set(gca, 'YDir','reverse')
caxis([-.4 .4])
c=colorbar;
title(c,gs03flma_velpt_variables{1,i}(3).units)
box on
title('GI03FLMA Meridional Velocities','fontweight','normal')

subplot(212)
hold on
for i = 1:4
    scatter((gs03flmb_velpt_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        (gs03flmb_velpt_variables{1,i}(9).data)*0.001,4,gs03flmb_velpt_variables{1,i}(3).data,'filled')
end
datetick('x',1,'keeplimits')
ylabel('dbar')
ylim([1600 2800])
set(gca, 'YDir','reverse')
caxis([-.4 .4])
c=colorbar;
title(c,gs03flmb_velpt_variables{1,i}(3).units)
box on
title('GI03FLMB Meridional Velocities','fontweight','normal')
