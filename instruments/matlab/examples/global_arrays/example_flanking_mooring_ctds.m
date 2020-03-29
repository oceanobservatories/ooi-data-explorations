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
%Specify metadata for Southern Ocean Flanking Mooring A
platform_name = 'GS03FLMA';
node = 'RISER';
instrument_class = 'CTD';
method = 'RecoveredInst';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);

for ii = 40:51
    disp(ii)
    if ii < 49
        dataset_name = strrep(uframe_dataset_name{1},'xx',num2str(ii));
        %Make M2M Call
        [nclist] = M2M_Call(dataset_name,start_date,end_date,options);
        %Get Data
        [gs03flma_ctd_variables{ii-39}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
    else
        dataset_name = strrep(uframe_dataset_name{2},'xx',num2str(ii));
        %Make M2M Call
        [nclist] = M2M_Call(dataset_name,start_date,end_date,options);
        %Get Data
        [gs03flma_ctd_variables{ii-39}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
    end
end

%%
%Specify metadata for Southern Ocean Flanking Mooring B
platform_name = 'GS03FLMB';
node = 'RISER';
instrument_class = 'CTD';
method = 'RecoveredInst';

%Get M2M URL
[uframe_dataset_name,variables] = M2M_URLs(platform_name,node,instrument_class,method);

for ii = 60:71
    disp(ii)
    if ii < 69
        dataset_name = strrep(uframe_dataset_name{1},'xx',num2str(ii));
        %Make M2M Call
        [nclist] = M2M_Call(dataset_name,start_date,end_date,options);
        %Get Data
        [gs03flmb_ctd_variables{ii-59}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
    else
        dataset_name = strrep(uframe_dataset_name{2},'xx',num2str(ii));
        %Make M2M Call
        [nclist] = M2M_Call(dataset_name,start_date,end_date,options);
        %Get Data
        [gs03flmb_ctd_variables{ii-59}, mtime, netcdfFilenames] = M2M_Data(variables, nclist);
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(211)
hold on
for i = 1:12
    plot((gs03flma_ctd_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),gs03flma_ctd_variables{1,i}(2).data,'k')
end
datetick('x',1)
ylabel(gs03flma_ctd_variables{1,i}(2).units)
ylim([2 12])
box on
title('GS03FLMA Seawater Temperature','fontweight','normal')

subplot(212)
hold on
for i = 1:12
    plot((gs03flmb_ctd_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),gs03flmb_ctd_variables{1,i}(2).data,'k')
end
datetick('x',1)
ylabel(gs03flmb_ctd_variables{1,i}(2).units)
ylim([2 12])
box on
title('GS03FLMB Seawater Temperature','fontweight','normal')


figure('units','normalized','outerposition',[0 0 1 1])
subplot(211)
hold on
for i = 1:12
    scatter((gs03flma_ctd_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        gs03flma_ctd_variables{1,i}(5).data,4,gs03flma_ctd_variables{1,i}(2).data,'filled')
end
datetick('x',1,'keeplimits')
ylabel(gs03flma_ctd_variables{1,i}(5).units)
ylim([0 1600])
set(gca, 'YDir','reverse')
caxis([0 12])
c=colorbar;
title(c,gs03flma_ctd_variables{1,i}(2).units)
box on
title('GS03FLMA Seawater Temperature','fontweight','normal')

subplot(212)
hold on
for i = 1:12
    scatter((gs03flmb_ctd_variables{1,i}(1).data)/60/60/24+datenum(1900,1,1),...
        gs03flmb_ctd_variables{1,i}(5).data,4,gs03flmb_ctd_variables{1,i}(2).data,'filled')
end
datetick('x',1,'keeplimits')
ylabel(gs03flmb_ctd_variables{1,i}(5).units)
ylim([0 1600])
set(gca, 'YDir','reverse')
caxis([0 12])
c=colorbar;
title(c,gs03flmb_ctd_variables{1,i}(2).units)
box on
title('GS03FLMB Seawater Temperature','fontweight','normal')
