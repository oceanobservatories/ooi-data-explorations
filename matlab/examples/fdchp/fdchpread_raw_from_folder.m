function [uw, U, sigH] = fdchpread_raw_from_folder(foldername)

    % Get file list
    filelist = dir(fullfile(foldername, '**/*.dat'));
    disp(["Processing ", num2str(length(filelist)), " files."]);
    parts = strsplit(foldername, 'uncabled/');
    parts = strsplit(parts{end}, '/');
    folder_end_name = parts{1};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read in the flux data file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strfind(folder_end_name,  'CP')
        lat=40.1334;           %Pioneer NES
        Rwaves=[-0.75 0 -5];   %
    elseif strfind(folder_end_name, 'GI')
        lat=59.9337;           %Irminger Sea
        Rwaves=[-0.75 0 -6];   %
    elseif strfind(folder_end_name,  'CE')
        lat=44.6393;           %Endurance
        Rwaves=[-0.5 -0.5 -5]; %
    else
        disp("Error: no settings available for folder ")
        return
    end
    L=12000; %Expected file length at 10 Hz
    bugfix=0;
    U=[];
    uw=[];
    sigH=[];
    
    for k=1:length(filelist)
        %****************************************
        % Read in raw data
        %****************************************
        if filelist(k).isdir
            continue
        end
        filename = [filelist(k).folder, '/', filelist(k).name]
        
        [td,dd] = fdchpread_raw(filename);
        rawdata=ones(21,L)*NaN;           %Preallocate
        rawdata=dd(1:21,1:L);
        
        %*****************************************
        % Compute flux data
        %*****************************************    
        [fluxes,Uearth,waveheight] = ProcessFDCHP_new(rawdata,L,bugfix,lat,Rwaves);
        uw=[uw -fluxes(1)];        %Fluxes: uw vw wT
        U=[U Uearth];              %Wind speed relative to earth  
        sigH=[sigH waveheight];    %Significant wave height
    end
    
    %******************************************
    % Plot it up
    %******************************************
    figure(1);clf
    plot(U,uw,'bo','markersize',5)
    xlabel('U (m/s)','fontsize',18)
    ylabel('-<uw> (m/s)','fontsize',18)
    axis([6.5 9 0.02 0.16])
    
    figure(2);clf
    plot(U,sigH,'bo','markersize',5)
    xlabel('U (m/s)','fontsize',18)
    ylabel('\sigma_H (m)','fontsize',18)
    axis([6.5 9 0 4])