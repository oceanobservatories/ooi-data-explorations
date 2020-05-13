function [variables, mtime, netcdfFilenames] = M2M_Data(variables, nclist, opendap)
%.. 2019-10-XX: CMRisien. original code.
%.. 2019-10-23: RADesiderio.
%..             the original version of this code was in a main program written
%..             by Craig Risien. I isolated it into this subroutine and used
%..             horzcat instead of specifying indices. 
%..       NOTE: (a) if the dataset is 1D, then the netcdf files contain
%..                 column vectors (time increases by increasing row index).
%..             (b) if the dataset is 2D, then the netcdf files contain 2D 
%..                 arrays in which time increases by increasing column index.
%.. 2019-11-04: RADesiderio. 
%..             (a) replaced nc_urls with nclist in the input argument list;
%..             (b) variable thredds_url transferred in from function M2M_Call;
%..             (c) added capability to read in the data from downloaded netcdf files
%..             (d) added calling argument opendap
%.. 2019-11-27: RADesiderio.
%..             (a) sorted nclist chronologically by deployment
%..             (b) appended deployment to variable list
%..             (c) read lat and lon from global attributes (non-mobile assets)
%
%.. INPUT: opendap is a switch. if it is 1 or true, the opendap call is made.
%..                             if it is 0 or false, nc files are downloaded.
%..                             DEFAULT IS TRUE.         
%
%
%.. USAGE:
%..         [variables, mtime, ~]   = M2M_Data(variables, nclist)   
%..         [variables, mtime, ~]   = M2M_Data(variables, nclist, true)   
%..
%..         [~, ~, netcdfFilenames] = M2M_Data(variables, nclist, false)
%
%.. NOTES:
%..    netcdf files containing data from surface moorings formerly had variables
%..    named 'lat' and 'lon'; now, they do not (except for some metbk and fdchp
%..    data streams). these netcdf files do have lat and lon mins and maxes 
%..    listed as global attributes (as do those metbk and fdchp exceptions).
%..    Therefore a section has been included to read in those global attribute
%..    values.
%
%.. this function was written to be run last in the sequence:
%..    M2M_URLs
%..    M2M_Call
%..    M2M_Data

if nargin==2
    opendap = true;
end

%.. additional variables to add, common to all instruments; enumerate here.
%.. .. lat and lon have been removed as netcdf variables for instruments
%.. .. not on mobile assets (gliders); will need to parse the netcdf
%.. .. files' global attributes separately to get lat and lon.
additionalVariables(1).name  = 'deployment';
additionalVariables(1).data  = [];
additionalVariables(1).units = 'unitless';
%additionalVariables(2).name  = 'fill in this space as appropriate';
%additionalVariables(2).data  = [];
%additionalVariables(2).units = 'fill in this space as appropriate';

if any(contains({variables.name}, 'lat'))
    isMobileAsset = true;
    %.. lat and lon are ordinary netcdf variables.
    %.. however, the mins and maxes are available from global attributes.
    %.. latMean and lonMean must be the LAST 2 VARIABLES.
    last2Variables(1).name  = 'latMean';
    last2Variables(1).data  = [];
    last2Variables(1).units = 'degrees_north';
    last2Variables(2).name  = 'lonMean';
    last2Variables(2).data  = [];
    last2Variables(2).units = 'degrees_east';
else
    isMobileAsset = false;
    %.. lat and lon must be the LAST 2 VARIABLES.
    last2Variables(1).name  = 'lat';
    last2Variables(1).data  = [];
    last2Variables(1).units = 'degrees_north';
    last2Variables(2).name  = 'lon';
    last2Variables(2).data  = [];
    last2Variables(2).units = 'degrees_east';
end

if isrow(variables)
    variables = [variables additionalVariables last2Variables];
elseif iscolumn(variables)
    variables = [variables; additionalVariables; last2Variables];
else
    error('Input structure ''variables'' must be a vector');
end

%.. determine filenames
%.. (a) for sorting nclist and ncfiles chronologically
%.. (b) in case ncfiles are downloaded (opendap==false)
netcdfFilenames = split(nclist, ["/", "\"], 2);
netcdfFilenames = netcdfFilenames(:, end);
%.. these begin with "deploymentXXXX_" (not 5X) where XXXX is deployment #
[netcdfFilenames, idx] = sort(netcdfFilenames);
%.. sort nclist so that data will be read in chronologically by deployment
nclist = nclist(idx);

%.. create a URL mapping of the files
thredds_url = "https://opendap.oceanobservatories.org/thredds/dodsC";
nc_urls = strcat(thredds_url, "/", nclist);  % string array

for i = 1:length(nclist)
    if opendap
        source = char(nc_urls(i));    % char required for R2018a (or use {i})
    else
        options = weboptions('Timeout',Inf);
        source = websave( ...
            netcdfFilenames(i), strrep(nc_urls{i,1}, 'dodsC', 'fileServer'),options);
        % char(source) *not* required for R2018a
    end
    %.. first variable is always time; pull it out of the following
    %.. loop so that nans can be assigned if a particular stream is
    %.. missing a variable (some nutnr cases)
    ncTime = ncread(source, variables(1).name);
    variables(1).data = horzcat(variables(1).data, ncTime');
    %.. if length(variables)==3 no error is thrown, and loop is not executed
    %.. (variables lat and lon are treated differently than the others);
    %.. subtract 2 because lat and lon are read from the global attributes
    for j = 2:length(variables)-2  
        try
            data = ncread(source, variables(j).name);
        catch
            %.. the following statement will cause an error at the first horzcat
            %.. statement if the missing variable is 2D
            data = nan(size(ncTime));
        end
        [~, col]=size(data);
        if col > 1
            %.. if the next line is referenced in a runtime error, see
            %.. comment immediately after the catch statement above.
            variables(j).data = horzcat(variables(j).data, data);
        else
            variables(j).data = horzcat(variables(j).data, data');
        end
    end
    %.. netdcf files for non-mobile assets have lat and lon information
    %.. under global attributes; these will be a little bit different
    %.. for each deployment.
    latMin = ncreadatt(source, '/', 'geospatial_lat_min');
    latMax = ncreadatt(source, '/', 'geospatial_lat_max');
    if ~isMobileAsset && latMin~=latMax
        disp('Warning! latMin and latMax don''t agree!');
    end
    lat = (latMin + latMax) / 2;
    data = zeros(size(ncTime)) + lat;
    variables(j+1).data = horzcat(variables(j+1).data, data');
    %
    lonMin = ncreadatt(source, '/', 'geospatial_lon_min');
    lonMax = ncreadatt(source, '/', 'geospatial_lon_max');
    if ~isMobileAsset && lonMin~=lonMax
        disp('Warning! lonMin and lonMax don''t agree!');
    end
    lon = (lonMin + lonMax) / 2;
    data = zeros(size(ncTime)) + lon;
    variables(j+2).data = horzcat(variables(j+2).data, data');
end
%.. first variable is always time; convert to matlab datenumber [day]
mtime=datenum(1900,1,1,0,0,0)+(variables(1).data/60/60/24);
