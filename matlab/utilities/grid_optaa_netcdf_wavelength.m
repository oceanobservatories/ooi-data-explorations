function grid_optaa_netcdf_wavelength(infile, outfile)
%.. 17-dec-2018 desiderio: initial code
%.. 11-jan-2019 desiderio: 
%..    collapsed the wavelength_a and wavelength_c dimensions into one 
%..    wavelength dimension.
%..    (1) this is possible because this program interpolates the 'a' and 'c'
%..        data onto a common wavelength grid.
%..    (2) for certain ERDDAP dataset types all variables must share a 
%..        common set of axes (dimensions), else some of the variables
%..        are not available on ERDDAP when the netcdf file is uploaded.
%.. 12-dec-2019 desiderio added tf_isnan construction at program end so that
%..                       interpolations won't fail when the OOINET netcdf
%..                       wavelength variables contain NaN values (!).
%.. 21-may-2020 desiderio 
%..             trapped out empty sets in the interp1 call for cases where
%..             all optaa vars are Nans (in these cases the CTD vars are
%..             not Nans).
%
%
%.. infile  = 'revamped' uframe OPTAA netcdf filename (nsif, mfn, or cspp);
%..           'revamped' means that the uframe netcdf file has been modified
%..                      by revamp_OOINet_optaa_netcdf_files.m.
%.. outfile = name of output netcdf file
%
%.. OPERATES ON PROCESSED OPTAA (WETLABS AC-S) UFRAME NETCDF FILES TO PUT 
%.. WAVELENGTH DATA ON A COMMON WAVELENGTH GRID.
%
%.. .. different factory dev files for the same ac-s instrument can have
%.. .. different numbers of wavelength channels; plus, the actual wavelength
%.. .. values change from factory cal to factory cal. In order to concatenate
%.. .. OPTAA netcdf data from consecutive deployments in ERDDAP a common
%.. .. wavelength grid is established by linear interpolation.

%.. NOTE
%.. Matlab 2018b does not support the NC_STRING netcdf datatype. If this is
%.. encountered this code will throw a runtime execution error. To date the
%.. uframe netcdf files encountered have not used this datatype and instead
%.. have used NC_CHAR.

gridded_wavelengths = 400:4:744;
nwvl_gridded = length(gridded_wavelengths);

schema = ncinfo(infile);

old_wavelength_dim_names = {'wavelength_a' 'wavelength_c'};
new_wavelength_dim_name  = {'wavelength'};

%.. remove one of the schema wavelength Dimensions substructures
tf = contains({schema.Dimensions.Name}, old_wavelength_dim_names{1});
schema.Dimensions(tf) = [];
%.. rename the other to 'wavelength'
tf = contains({schema.Dimensions.Name}, old_wavelength_dim_names{2});
schema.Dimensions(tf).Name = new_wavelength_dim_name{1};
%.. and change the length of the wavelength dimension to 
%.. the number of gridded wavelengths
[schema.Dimensions(tf).Length] = nwvl_gridded;

%.. same with the wavelength_a and wavelength_c variables ...
idx_wvl_var = find(contains({schema.Variables.Name}, old_wavelength_dim_names));
%.. delete the substructure for the higher indexed wavelength_x variable;
%.. therefore the index for the lower one will not change 
schema.Variables(idx_wvl_var(2)) = [];
%.. keep the other:
%.. change the variable name 
schema.Variables(idx_wvl_var(1)).Name = new_wavelength_dim_name{1};
%.. change the name of the variable's dimension
schema.Variables(idx_wvl_var(1)).Dimensions.Name = new_wavelength_dim_name{1};
%.. change the length value for the variable's dimension
schema.Variables(idx_wvl_var(1)).Dimensions.Length = nwvl_gridded;
%.. reset the Size and ChunkSize
schema.Variables(idx_wvl_var(1)).Size      = nwvl_gridded;
schema.Variables(idx_wvl_var(1)).ChunkSize = nwvl_gridded;
%.. delete the comment and long_name Attributes which are now obsolete
tf_comment = contains({schema.Variables(idx_wvl_var(1)).Attributes.Name}, ...
    {'comment' 'long_name'});
schema.Variables(idx_wvl_var(1)).Attributes(tf_comment) = [];

%.. get the indices of the 1D and 2D variables
idx_var1D = find(cellfun('length', {schema.Variables.Dimensions})==1);
idx_var2D = find(cellfun('length', {schema.Variables.Dimensions})==2);

%.. reset the dimension Names and Lengths of the 2D variables, and, 
%.. .. their Sizes and ChunkSizes
for ii = 1:length(idx_var2D)
    tf = contains({schema.Variables(idx_var2D(ii)).Dimensions.Name}, ...
        old_wavelength_dim_names);
    %.. dimension name
    schema.Variables(idx_var2D(ii)).Dimensions(tf).Name = ...
        new_wavelength_dim_name{1};
    %.. dimension length
    schema.Variables(idx_var2D(ii)).Dimensions(tf).Length = nwvl_gridded;
    %.. size
    schema.Variables(idx_var2D(ii)).Size(tf) = nwvl_gridded;
    %.. chunksize
    schema.Variables(idx_var2D(ii)).ChunkSize(tf) = nwvl_gridded;
end

%.. SET UP NETCDF OUTFILE
if ~strcmp(outfile(end-2:end), '.nc')
    outfile = [outfile '.nc'];
end
%.. MAKE SURE THAT THIS FILE DOES NOT ALREADY EXIST!
%.. by design the ncwriteschema function *appends* schema to the file,
%.. if it exists; in this code we want a completely new file.
if isfile(outfile)
    action2take = 'Either remove or rename it, or use a different outfilename.';
    disp(' ');
    error('Outfilename ''%s'' exists as ''%s''\n%s\n\n\n', ...
        outfile, which(outfile), action2take);
end
%.. create empty netcdf file with new schema.
ncwriteschema(outfile, schema);

%.. write out the gridded wavelength values
ncwrite(outfile, new_wavelength_dim_name{1}, gridded_wavelengths');

%.. write out 1D data variables excluding wavelength data;
%.. schema no longer has the old wavelength variable names, but it 
%.. does have the new name.
varname = setdiff({schema.Variables(idx_var1D).Name}, ...
    new_wavelength_dim_name);  %#ok
for ii = 1:length(varname) 
    tmp = ncread(infile, varname{ii});
    %.. in case this code is adapted to use on raw uframe netcdf files,
    %.. trap out infinities and replace with nans.
    tmp(isinf(tmp(:))) = nan;
    ncwrite(outfile, varname{ii}, tmp);
end

%.. read in actual acs wavelengths 
wvl_a = ncread(infile, 'wavelength_a');
wvl_c = ncread(infile, 'wavelength_c');
%.. some OOINET datasets use NaN values to pad the wavelength record!
%.. eliminate these so that interpolations will run.
tf_isnan = isnan(wvl_a);
wvl_a(tf_isnan) = [];
wvl_c(tf_isnan) = [];

%.. now read in and grid each 2D data variable using linear interpolation
varname = {schema.Variables(idx_var2D).Name};
for ii = 1:length(varname)
    %.. these 2D arrays are sized as (nwvl_acs x time) 
    tmp = ncread(infile, varname{ii});
    %.. eliminate bogus rows if there were Nan values in the wavelength record
    tmp(tf_isnan, :) = [];
    %.. in case this code is adapted to use on raw uframe netcdf files,
    %.. trap out infinities and replace with nans.
    tmp(isinf(tmp(:))) = nan;
    %.. determine whether 'a' or 'c' wavelengths should be used
    if contains(varname{ii}, {'a_', 'abs'})
        wvl = wvl_a;
    elseif contains(varname{ii}, {'c_', 'att'})
        wvl = wvl_c;
    else
        error('Could not assign ''a'' or ''c'' to a 2D spectral variable.');
    end
    %.. interpolate
    %.. .. trap out empty sets
    if isempty(wvl), break, end
    if isempty(tmp), continue, end
    tmp = interp1(wvl, tmp, gridded_wavelengths');
    %.. write out
    ncwrite(outfile, varname{ii}, tmp);
end
