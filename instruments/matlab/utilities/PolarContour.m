function [px, py] = PolarContour(dirSpec, freq)
%PolarContour creates a polar/contour plot for the presentation of
%directional data. This Plot was orginally designed to plot directional 
%frequency spectra from ocean waves.
%   dirSpec - A matrix of directional frequency data, with frequency data
%             propagating along the x-axis (horizontal) and direction across
%             the y-axis (vertical).
%   freq    - An index of frequency intervals in a vertical list.
%

%Open a new figure
fig=figure;                                                                 

%Defines directions bins based on dirSpec x dimentions (in radians) 
temp=360-(360/size(dirSpec,1));
a=linspace(0,temp,size(dirSpec,1));                                          
dirs=degtorad(a);
 
%Creates artificial dataset to plot in polar cordinates
[df,ddir]=meshgrid(freq,a); 
ddir=degtorad(ddir);                                                        
[px,py]=pol2cart(ddir,df);                                                  
h = polar(px,py);                                                           

%Removes artificial dataset but keeps axis
delete(h);                                                                  
hold on                                                                     

%test the code below
%px=[px;px(1,:)];
%py=[py;py(1,:)];
%dirSpec=[dirSpec;dirSpec(1,:)];

%Plots contour graph on to polar axis
contour(px,py,dirSpec,20);                                                       

%Plot aesthetics
%Roates and flips plot to enable 0 degrees = North
view([90 -90])
%Sets axis ticks for diffined values (frequency values in this case)
set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5])                                      
%colorbar('vert');                                                           
%ylabel('Direction [degrees] / Frequency [Hz]');                             
%xlabel('m^2s / deg');

end

