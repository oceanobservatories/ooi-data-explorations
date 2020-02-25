function [variables] = M2M_SPKIR(netcdfFilenames)

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

%.. read in downwelling data: get variable names (vars)
nc_struct = ncinfo(processed_file{1});
vars = {nc_struct.Variables.Name}';
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
% make sure bounds outside actual data ranges
fractionalDay_T000730 = 7/60/24 + 30/60/60/24;
bin_first = floor(min(mtime-1)) + fractionalDay_T000730;
bin_last  = floor(max(mtime+1)) + fractionalDay_T000730;
bin_edges = (bin_first*24*60:15:bin_last*24*60)/24/60; %15 min bins
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
