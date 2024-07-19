clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the flux data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PA=0;EA=0;IS=1;
if PA
    lat=40.1334;           %Pioneer NES
    Rwaves=[-0.75 0 -5];   %
elseif IS
    lat=59.9337;           %Irminger Sea
    Rwaves=[-0.75 0 -6];   %
else
    lat=44.6393;           %Endurance
    Rwaves=[-0.5 -0.5 -5]; %
end
L=12000;
bugfix=0;
U=[];
uw=[];
sigH=[];

for incr=1:23
    %****************************************
    % Read in raw data
    %****************************************
    
    if incr<=9
        filename = ['fdchp_20200818_0' num2str(incr) '0200.dat'];
    else
        filename = ['fdchp_20200818_' num2str(incr) '0200.dat'];
    end
    [td,dd] = fdchpread_raw(filename);
    L=12000;                          %Expected length at 10 Hz
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

