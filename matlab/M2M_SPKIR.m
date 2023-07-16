function [variables] = M2M_SPKIR(netcdfFilenames, tf_medianFilter)
%.. 2019-12-18: CMRisien. original code
%.. 2020-06-15: RADesiderio. 
%..             (a) added optional calling argument tf_medianFilter; if omitted,
%..                 tf_medianFilter is set to true to preserve original behavior.
%..             (b) added calculation of binsize from time data.
%.. 2020-06-19: RADesiderio. include pressure data on output if it exists.

if nargin==1, tf_medianFilter=true; end

netcdfFilenames = string(netcdfFilenames);  % so that character vector will work
ooifile = netcdfFilenames(contains(netcdfFilenames, 'SPKIR'));
processed_file = strcat('processed_', ooifile);
for i = 1:length(ooifile)
    revamp_OOINet_spkir_netcdf_files(ooifile{i},processed_file{i});
end

mtime=[];
for i = 1:length(processed_file)
    time=ncread(processed_file{i},'time')/60/60/24;
    time=time+datenum(1900,1,1,0,0,0);
    mtime = horzcat(mtime, time');
end
variables.mtime = mtime;

%.. get the SPKIR wavelengths
variables.wavelength = ncread(processed_file{1}, 'wavelength');

%.. get variable names (vars)
nc_struct = ncinfo(processed_file{1});
vars = {nc_struct.Variables.Name}';

%.. if pressure data are present, 
%.. add the data to the structure as a column vector
if any(strcmp(vars, 'pressure'))
    pvarName = 'pressure';
elseif any(strcmp(vars, 'pressure_depth'))
    pvarName = 'pressure_depth';
elseif any(strcmp(vars, 'int_ctd_pressure'))
    pvarName = 'int_ctd_pressure';
else
    pvarName = '';
end
if ~isempty(pvarName)
    pressure=[];
    for i = 1:length(processed_file)
        ppp=ncread(processed_file{i}, pvarName);
        pressure = vertcat(pressure, ppp(:));  %#ok
    end
    variables.pressure = pressure;
end

%.. read in downwelling data: 
vars = vars(contains(vars, 'Ed_'));
%Cat the data
for ii = 1:length(vars)
    varII=[];
    for jj = 1:length(processed_file)
        tmpJJ=ncread(processed_file{jj}, vars{ii});
        varII = horzcat(varII, tmpJJ');
    end
    variables.(vars{ii}) = varII;
end

%%
%Bin raw data and apply median
%.. for profiler data set calling argument tf_medianFilter to false
if tf_medianFilter
    %.. determine binsize in minutes
    dTime = 60 * 24 * diff(mtime);
    %.. remove diffs during data acquisition within bursts
    dTime(dTime<1) = [];  % daq rates less than 1 minute per point   
    %.. round the most common diff up to the nearest 15 minutes
    binsize = 15 * ceil(nanmedian(dTime)/15);
    disp(['M2M_SPKIR: binsize calculated to be ' num2str(binsize) ' minutes.']);
    variables.binsizeMinutes = binsize;
    % make sure bounds outside actual data ranges
    fractionalDay_halfBinsize = 0.5 * binsize/60/24;
    bin_first = floor(min(mtime)) - fractionalDay_halfBinsize;
    bin_last  = ceil(max(mtime))  + fractionalDay_halfBinsize;
    bin_edges = (bin_first*24*60:binsize:bin_last*24*60)/24/60;
    %
    ind = discretize(mtime,bin_edges);
    % time
    mtime_binned = accumarray(ind',mtime,[],@nanmedian);
    mtime_binned(mtime_binned==0)=nan;
    variables.mtime_binned = mtime_binned';
    
    % binned wavelength channels
    binned_vars = strcat(vars, '_binned');
    for ii = 1:length(vars)
        binnedII = accumarray(ind',variables.(vars{ii}),[],@nanmedian);
        binnedII(binnedII==0)=nan;
        variables.(binned_vars{ii}) = binnedII';
    end
end
