function schema = revamp_OOINet_optaa_netcdf_files(infilename, outfilename, blankingTimeSec, tf_addSpectral1D)
%.. 13-dec-2018 desiderio Matlab 2018b initial code
%.. 12-dec-2019 desiderio
%..             (a) added tf_isnan construction at program end so that
%..                 interpolations won't fail when the OOINET netcdf
%..                 wavelength variables contain NaN values (!).
%..             (b) exclude warm-up time in burst median calculation.
%.. 13-dec-2019 desiderio
%..             added switches to control median filtering and spectral1D variables
%..             (a) burstMedian filtering can be disabled, or enabled with a 
%..                 default or user-supplied blanking time.
%..             (b) a switch has been added to control whether spectral 1D time
%..                 series (akin to ac-9 data) are to be added.
%.. 17-dec-2019 desiderio added 660nm to spectral1D wavelength set
%.. 22-may-2020 desiderio
%..             BEP netcdf varname seawater_temperature changed to temperature
%.. 27-may-2020 desiderio
%..             (a) trapped out case where blanking time >= burst time would result
%..                 in index out of bounds error.
%..             (b) fixed mechanism of singleton dimension elimination after median
%..                 filtering for case when nbin = 1: squeeze.m results in column
%..                 vector instead of rowvector when size(array) = [1 1 N] so it
%..                 was replaced with filteredArray assignment statements
%.. 29-may-2020 desiderio
%..             formerly when the number of points in a burst was less than the number 
%..             of blanking points a time value of nan resulted. now replaces these 
%..             values of nan with times derived from unblanked data.
%
%.. USAGE (DEFAULTS)
%
%..     This function can be called with 2, 3, or 4 calling arguments.   
%..     nargin == 2:
%..          BURST MEDIAN FILTERING IS ENABLED USING THE DEFAULT BLANKING TIME.
%..          tf_addSpectral1D is set to false.
%..     nargin == 3:
%..          burst median filtering is controlled by the value of the 3rd argument.
%..          tf_addSpectral1D is set to false.
%..     nargin == 4:
%..          burst median filtering is controlled by the value of the 3rd argument.
%..          addition of 1D time series is controlled by the value of the 4th argument.
%
%
%.. CALLING ARGUMENTS
%
%.. infile  = OOINet (from uframe) netcdf filename containing surface mooring
%..           or cspp OPTAA (Seabird\WETLabs\ac-s) data.
%
%.. outfile = name of output netcdf file containing revamped netdcf variables,
%..           dimensions, and data
%
%.. blankingTimeSec serves multiple purposes:
%..         = nan, then data is not burst median filtered (see Notes)
%..         = [],  then burst median filtering is enabled using a default value
%..                set by the variable blankingTimeSecDefault at code top for
%..                the blanking time (during which the instrument is warming up
%..                so that the corresponding data are excluded from the filtering).
%..                This default value has been set at 60 seconds (12/13/2019).
%..         = x,   where x >=0, then burst median filtering is enabled using 
%..                the input value x as the blanking time in seconds. A value 
%..                of 0 result in no blanking and therefore no exclusion of 
%..                data - all data are filtered.
%
%.. tf_addSpectral1D - controls whether ac-9-like time series are added to
%..                    the data. Originally added to provide discrete (1D) 
%..                    wavelength channels for absorption and beam attenuation
%..                    time series data for visualization in Highcharts.
%..             The wavelengths selected are given in the variable wavelengths_1D. 
%..         = true  or 1: additional times series variables added.           
%..         = false or 0: no action is taken (default).
%
%
%.. ACTIONS
%
%.. reprocesses OPTAA (WETLabs ac-s) UFrame netcdf files
%
%.. (0) finds the number of ac-s wavelength channels
%.. (1) deletes extraneous dimensions
%.. (2) deletes extraneous variables
%.. (3) declares variable wavelength_a as a dimension (-> coordinate variable)
%.. (4) declares variable wavelength_c as a dimension (-> coordinate variable)
%.. (5) changes dimension of variable wavelength_a to wavelength_a
%.. (6) changes dimension of variable wavelength_c to wavelength_c
%.. (7) changes wavelength dimension of 2D absorption variables to wavelength_a 
%.. (8) changes wavelength dimension of 2D beam c variables to wavelength_c 
%.. (9) replaces 'obs' dimension with 'time' to make time a coordinate variable
%.. (A) renames some variables to more standardized names
%.. (B) renames above renamed variables in ancillary_variables attribute of 
%..     beam_attenuation and optical_absorption variables
%.. (C) writes units of temperature as 'degree_C'
%.. (D) changes fill values of L2 data products to NaN.
%.. (E) checks and warns if 'UNSUPPORTED DATATYPE' is found in matlab's schema
%.. (F) optionally median filters the data in each burst
%.. (G) optionally splits out ~11 wavelength channels for optical_absorption
%.. (H) optionally splits out ~11 wavelength channels for beam_attenuation
%.. (I) writes out 1 netcdf file using the new schema 

%.. NOTES
%
%****************************************************************************
%.. DO NOT PROCESS CSPP OPTAA DATA USING BURST MEDIAN FILTERING. 
%.. Deactivate this option by calling the function with nan as the 3rd argument.
%
%.. If burst medain filtering is used on cspp data then each profile will have
%.. all of its OPTAA data channels median filtered to one value each per profile.
%****************************************************************************
%
%.. Burst Median Filtering
%..
%.. OPTAAs installed on surface moorings are scheduled to cycle power on and
%.. off throughout the duration of the deployment. Typically data are acquired
%.. at 4 Hz for about 3 minutes then the instrument is powered off for the next
%.. 12 minutes or longer, depending on the platform and power availability. A
%.. blanking period is often implemented to exclude data acquired just after
%.. power on when the instrument is warming up.
%..
%.. The burst median filter algorithm requires a time parameter to determine 
%.. burst groupings. This value must be greater than the time interval between
%.. data points within a burst and less than the minimum time between the last
%.. datum of a burst and the first datum of the next burst. For some OPTAA 
%.. deployments on surface moorings, the latter value is about 15-3=12 minutes.
%.. For general applicability this value is set to 300 seconds at code top as: 
%..      minSleepTimeBetweenBurstsSec = 300;
%
%
%.. Nan values in netcdf wavelength variables 
%..
%.. Certain machine to machine time range requests for OOINet data will result 
%.. in nan values in the (input) netcdf wavelength vectors. For this reason the
%.. tf_isnan construction was implemented at the end of the code. The presence
%.. of nans in the wavelength variables in the output file will prevent plots
%.. using wavelength as an axis when the outfile is opened in panoply.
%..
%.. The next program in the processing sequence is grid_optaa_netcdf_wavelength.m.
%.. On output this code automatically eliminates nans in the wavelength records
%.. because the 2D spectral data is interpolated onto a set wavelength grid.
%
%
%.. Matlab 2018b and NetCDF compatability
%..
%.. Matlab 2018b does not support the NC_STRING netcdf datatype. If this is
%.. encountered this code will throw a runtime execution error. To date the
%.. uframe netcdf files encountered have not used this datatype and instead
%.. have used NC_CHAR.

minSleepTimeBetweenBurstsSec = 300;
%************************************************************************
%***** PROCESS INPUT ARGUMENT LIST **************************************
%************************************************************************
blankingTimeSecDefault = 60;
if nargin==2
    tf_addSpectral1D = false;
    tf_burstMedian   = true;
    blankingTimeSec  = blankingTimeSecDefault;
elseif nargin==3 || nargin==4
    if isempty(blankingTimeSec)
        tf_burstMedian   = true;
        blankingTimeSec  = blankingTimeSecDefault;
    elseif isnan(blankingTimeSec)
        tf_burstMedian   = false;
    elseif ~isnumeric(blankingTimeSec) || blankingTimeSec < 0
        error('Blanking time must be numeric (and non-negative).')
    else
        tf_burstMedian   = true;
    end
    
    if nargin==3
        tf_addSpectral1D = false;
    end
    
    if tf_addSpectral1D  
        disp('1D spectral time series to be added.');
        wavelengths_1D = [412 443 490 520 555 585 620 650 660 676 715];
        nwvl_1D = length(wavelengths_1D);
    end
end


%************************************************************************
%***** CONDITION DIMENSIONS AND VARIABLES FROM INPUT NETCDF FILE ********
%************************************************************************

%.. old dimension name(s) to keep from infile
dims2keep = {'obs'};
%.. this will be renamed as:
old_abscissa_dim_name = 'obs';
new_abscissa_dim_name = 'time';

%.. variables (uframe names) to keep from infile. 
%.. .. some of these are unique to surface moorings, some to cspps.
%** NOTE AGAIN THAT THIS CODE WILL CONSIDER EACH CSPP PROFILE TO BE ONE
%** BIN OF DATA AND MEDIAN FILTER THE OPTAA DATA ACCORDINGLY.
vars2keep = {
    'a_reference_counts'
    'a_reference_dark_counts'
    'a_signal_counts'
    'a_signal_dark_counts'
    'beam_attenuation'
    'c_reference_counts'
    'c_reference_dark_counts'
    'c_signal_counts'
    'c_signal_dark_counts'
    'elapsed_run_time'
    'external_temp_raw'
    'internal_temp_raw'
    'int_ctd_pressure'
    'lat'
    'lon'
    'on_seconds'
    'optical_absorption'
    'practical_salinity'
    'pressure_depth'
    'salinity'
    'seawater_temperature'
    'temp'
    'temperature'
    'time'
    'wavelength_a'
    'wavelength_c'
    };
%.. variables to be renamed in outfile
oldVariableName = {
    'external_temp_raw'     % all
    'internal_temp_raw'     % all
    'on_seconds'            % cspp name
    'practical_salinity'    % bep and surface mooring name
    'seawater_temperature'  % bep name
    'temp'                  % surface mooring name
    };
newVariableName = {
    'external_acs_temperature_counts' 
    'internal_acs_temperature_counts'
    'elapsed_run_time'      % bep and surface mooring name
    'salinity'              % cspp name
    'temperature'           % cspp name        
    'temperature'           % (twice to change both bep and surface mooring)        
    };

schema = ncinfo(infilename);

%.. find number of acs wavelengths
idx = find(strcmp({schema.Dimensions.Name}, 'wavelength'));
nwvl_acs = schema.Dimensions(idx).Length;
%.. and save substructure for writing in new dimensions
dim_wvl = schema.Dimensions(idx);

%.. delete extraneous dimensions
alldims = {schema.Dimensions.Name};
schema.Dimensions(~ismember(alldims, dims2keep)) = [];

%.. delete the variables not listed in vars2keep
allvars = {schema.Variables.Name};
schema.Variables(~ismember(allvars, vars2keep)) = [];

%.. add new dimensions
dim_wvl.Name = 'wavelength_a';
schema.Dimensions(end+1) = dim_wvl;
dim_wvl.Name = 'wavelength_c';
schema.Dimensions(end+1) = dim_wvl;

%.. and change the dimensions of the corresponding 1D variable of the same name
tf = strcmp({schema.Variables.Name}, 'wavelength_a');
schema.Variables(tf).Dimensions.Name = 'wavelength_a';
tf = strcmp({schema.Variables.Name}, 'wavelength_c');
schema.Variables(tf).Dimensions.Name = 'wavelength_c';

%.. find the indices of the 1D and 2D variables
ndims_of_vars = cellfun('length', {schema.Variables.Dimensions});
idx_var1D = find(ndims_of_vars==1);  % will be used in 1D median filtering
idx_var2D = find(ndims_of_vars==2);
%.. change the dimension declarations for the corresponding 2D variables
for ii = idx_var2D
    %.. need to change 'wavelength' to either 'wavelength_a' or 'wavelength_c'
    idx_wavelength = find(strcmp({schema.Variables(ii).Dimensions.Name}, ...
        'wavelength'));
    if strcmp(schema.Variables(ii).Name(1:2), 'a_')
        schema.Variables(ii).Dimensions(idx_wavelength).Name = 'wavelength_a';
    elseif strcmp(schema.Variables(ii).Name(9:11), 'abs')
        schema.Variables(ii).Dimensions(idx_wavelength).Name = 'wavelength_a';
    elseif strcmp(schema.Variables(ii).Name(1:2), 'c_')
        schema.Variables(ii).Dimensions(idx_wavelength).Name = 'wavelength_c';
    elseif strcmp(schema.Variables(ii).Name(1:8), 'beam_att')
        schema.Variables(ii).Dimensions(idx_wavelength).Name = 'wavelength_c';
    end
end

%.. change the abscissa dimension name(s) from old to new
tf = strcmp({schema.Dimensions.Name}, old_abscissa_dim_name);
schema.Dimensions(tf).Name = new_abscissa_dim_name;

%.. change the variables' abscissa dimension assignations from old to new
for ii = 1:length(schema.Variables)
    for jj = 1:length(schema.Variables(ii).Dimensions)
        if strcmp(schema.Variables(ii).Dimensions(jj).Name, ...
                old_abscissa_dim_name)
            schema.Variables(ii).Dimensions(jj).Name = new_abscissa_dim_name;
        end
    end
end

%.. change the listed 'oldVariableName's to 'newVariableName's
%.. .. these lists are updated according to the type of OPTAA file being 
%.. .. processed by the statements in the if\else construction so that
%.. .. the old data can be read in from the old file and written out
%.. .. to the new file.
%..
%.. .. Because of the deletions in the 'else' section, iterate the forloop
%.. .. in the reverse direction. 
for ii = length(oldVariableName):-1:1 
    tf = strcmp({schema.Variables.Name}, oldVariableName{ii});
    if any(tf) 
        %.. if oldname found, change it 
        schema.Variables(tf).Name = newVariableName{ii};
    else
        %.. else, delete it and its replacement name;
        %.. this construction enables both cspp and surface mooring
        %.. optaa netcdf files to be run with this code
        oldVariableName(ii) = [];
        newVariableName(ii) = [];
    end
end

%.. some of the variables whose names have been changed show up in
%.. ancillary variable attributes; change these to new names
for ii = 1:length(schema.Variables)
    idx = find(strcmp({schema.Variables(ii).Attributes.Name}, ...
        'ancillary_variables'));
    if ~isempty(idx)
        str = schema.Variables(ii).Attributes(idx).Value;
        %.. because 'temp' appears in several variable names, isolate by
        %.. adding a comma to the start and end of the list, then search
        %.. and replace ',temp,'
        str = strtrim(str);  % make sure there are no spaces.
        str = [',' str ','];  %#ok
        for jj = 1:length(oldVariableName)
            str = strrep(str, ...
                [',' oldVariableName{jj} ','], [',' newVariableName{jj} ',']);
        end
        str(1)   = [];
        str(end) = [];
        schema.Variables(ii).Attributes(idx).Value = str;
    end
end

%.. find units for 'temperature' which has been read in to matlab as
%.. 'UNSUPPORTED DATATYPE' because Uframe used the degree symbol, 
%.. which matlab does not support. 
tf_name = strcmp({schema.Variables.Name}, 'temperature');
tf_units = strcmp({schema.Variables(tf_name).Attributes.Name}, 'units');
%.. can't use 'C', which stands for Coulomb. 
schema.Variables(tf_name).Attributes(tf_units).Value = 'degree_C';

%.. change the fill value attribute of the L2 data products to NaN.
tf_name = ismember({schema.Variables.Name}, ...
    {'optical_absorption' 'beam_attenuation'});
[schema.Variables(tf_name).FillValue] = deal(NaN);

%.. send out a warning if 'UNSUPPORTED DATATYPE' shows up in either
%.. the Variables datatype entries or the Variables.Attributes units
%.. entries.
datatypes = extractfield(schema.Variables, 'Datatype');
tf_unsupported = strcmp(datatypes, 'UNSUPPORTED DATATYPE');
if any(tf_unsupported)
    disp('Warning:');
    disp(' ''UNSUPPORTED DATATYPE'' found in Variables Datatype field:');
    disp({schema.Variables(tf_unsupported).Name}');
end
%.. now look in the Attributes Value fields (which includes values for units)
for ii = 1:length(schema.Variables)
    values = {schema.Variables(ii).Attributes.Value};
    tf_unsupported = strcmp(values, 'UNSUPPORTED DATATYPE');
    if any(tf_unsupported)
        disp('Warning:');
        disp(' ''UNSUPPORTED DATATYPE'' found in Variables Attributes field.');
        disp(['Variable Name: ' schema.Variables(ii).Name]);
        name  = schema.Variables(ii).Attributes(tf_unsupported).Name;
        disp(['Attribute    : ' name]);    
    end
end

%****************************************************************************
if tf_addSpectral1D
    %************************************************************************
    %.. ADD NETCDF INFO FOR 1D VARIABLES TO SCHEMA
    %************************************************************************
    
    %.. construct spectral 1D variable names
    abs_1D_varname(1:nwvl_1D)   = {''};
    beam_c_1D_varname(1:nwvl_1D) = {''};
    for ii = 1:nwvl_1D
        wave_char = [num2str(wavelengths_1D(ii)) 'nm'];
        abs_1D_varname{ii}   = ['abs_' wave_char];
        beam_c_1D_varname{ii} = ['beam_c_' wave_char];
    end
    
    %.. add netcdf schema specifications for 1D absorption variables
    %.. .. isolate corresponding 2D Variables structure to serve as a template
    tf = strcmp({schema.Variables.Name}, 'optical_absorption');
    abs_template = schema.Variables(tf);
    %.. reset dimensions from 2D to 1D by eliminating wavelength dimension
    tf = strcmp({abs_template.Dimensions.Name}, 'wavelength_a');
    abs_template.Dimensions(tf) = [];
    %.. reset Size 2D -> 1D
    abs_template.Size(abs_template.Size==nwvl_acs) = [];
    %.. reset ChunkSize 2D -> 1D
    abs_template.ChunkSize(abs_template.ChunkSize==nwvl_acs) = [];
    %.. add the nwvl variables to the schema structure
    for ii = 1:nwvl_1D
        abs_template.Name = abs_1D_varname{ii};
        schema.Variables(end+1) = abs_template;
    end
    
    %.. add netcdf schema specifications for 1D beam attenuation variables
    %.. .. isolate corresponding 2D Variables structure to serve as a template
    tf = strcmp({schema.Variables.Name}, 'beam_attenuation');
    beam_c_template = schema.Variables(tf);
    %.. reset dimensions from 2D to 1D by eliminating wavelength dimension
    tf = strcmp({beam_c_template.Dimensions.Name}, 'wavelength_c');
    beam_c_template.Dimensions(tf) = [];
    %.. reset Size 2D -> 1D
    beam_c_template.Size(beam_c_template.Size==nwvl_acs) = [];
    %.. reset ChunkSize 2D -> 1D
    beam_c_template.ChunkSize(beam_c_template.ChunkSize==nwvl_acs) = [];
    %.. add the nwvl variables to the schema structure
    for ii = 1:nwvl_1D
        beam_c_template.Name = beam_c_1D_varname{ii};
        schema.Variables(end+1) = beam_c_template;
    end
    %************************************************************************
    %.. END:  ADD NETCDF INFO FOR 1D VARIABLES TO SCHEMA
    %************************************************************************
end  % tf_addSpectral1D
%****************************************************************************

%************************************************************************
%***** SET UP NETCDF OUTFILE ********************************************
%************************************************************************
if ~strcmp(outfilename(end-2:end), '.nc')
    outfilename = [outfilename '.nc'];
end
%.. MAKE SURE THAT THIS FILE DOES NOT ALREADY EXIST!
%.. by design the ncwriteschema function *appends* schema to the file, if it
%.. exists; in this code we want a completely new file.
if isfile(outfilename)
    action2take = 'Either remove or rename it, or use a different outfilename.';
    disp(' ');
    error('Outfilename ''%s'' exists as ''%s''\n%s\n\n\n', ...
        outfilename, which(outfilename), action2take);
end

%.. create empty netcdf file with new schema.
ncwriteschema(outfilename, schema);
%****************************************************************************

%************************************************************************
%***** WRITE OUT BASE DATA **********************************************
%************************************************************************
if ~tf_burstMedian    % ........... NO FILTERING ........................
    %.. write out (kept) old data first,
    %.. but exclude newly created spectral 1D variable names if they exist
    if tf_addSpectral1D
        varnames = ...
            setdiff({schema.Variables.Name}, [abs_1D_varname beam_c_1D_varname]);
    else
        varnames = {schema.Variables.Name};
    end
    for ii = 1:length(varnames)
        %.. if necessary change new varnames back to old to read data from infile
        outfile_varname = varnames{ii};
        tf = ismember(newVariableName, outfile_varname);
        if any(tf)
            infile_varname = oldVariableName{tf};
        else
            infile_varname = outfile_varname;
        end
        tmp = ncread(infilename, infile_varname);
        %.. trap out infinities and replace with nans.
        tmp(isinf(tmp(:))) = nan;
        ncwrite(outfilename, outfile_varname, tmp);
    end
else    % .............. BURST MEDIAN FILTERING ........................
    %.. write out variables that are not to be median filtered
    var_doNotFilter = {'wavelength_a' 'wavelength_c'};
    for ii = 1:length(var_doNotFilter)
        tmp = ncread(infilename, var_doNotFilter{ii});
        ncwrite(outfilename, var_doNotFilter{ii}, tmp);
    end
    %.. MEDIAN FILTER THE DATA IN EACH BURST
    time = ncread(infilename, 'time');
    npts = length(time);
    %.. set up burst endpoints and bin sizes.
    %.. append a time interval large enough to set the last endpoint
    dt = [diff(time); 2*minSleepTimeBetweenBurstsSec];
    %.. prepend 0 to set starting 'endpoint'
    endpoint = [0; find(dt > minSleepTimeBetweenBurstsSec)];
    burstsize = diff(endpoint);
    binsize = max(burstsize);
    nbins = length(burstsize);
    %.. set up 'blanking' - how many points to skip till the acs is warmed up.
    timeBetweenPoints = nanmedian(dt);  % [seconds]
    nptsToSkip        = floor(blankingTimeSec/timeBetweenPoints) + 1;
    %.. nptsToSkip=0 is OK, 1:0=[];
    %.. as an array index this stride results in no action.
    
    %.. MEDIAN FILTER 1D VARIABLES
    %.. .. but not spectral 1D variables = f(wavelength) if these are to be created;
    %.. .. if so, they will be created after interpolation of the 2D data.
    %.. indices in schema of 1D variables to be filtered:
    idx_var_not_to_filter = ...
        find(ismember({schema.Variables.Name}, var_doNotFilter));
    idx_var1D_to_filter = setdiff(idx_var1D, idx_var_not_to_filter);
    nvars = length(idx_var1D_to_filter);
    
    %.. track the names of the 1D filtered variables
    var1D_filtered_new_names(1:nvars) = {''};
    %.. read in 1D data from infile and place in a data array
    data_array = NaN(npts, nvars);
    for ii = 1:nvars
        varname = schema.Variables(idx_var1D_to_filter(ii)).Name;
        %.. save it here to make sure new names are saved
        var1D_filtered_new_names{ii} = varname;  
        %.. if necessary change new varnames back to old to read in data
        tf = ismember(newVariableName, varname);  % tf is a logical vector
        if any(tf)
            varname = oldVariableName{tf};
        end
        data_array(:, ii) = ncread(infilename, varname);
    end
    
    %.. change infinities, which are often not compatible with netcdf tools, to nan
    data_array(isinf(data_array(:))) = nan;
    %.. set up for fast median filtering
    arr2filter = NaN(binsize*nbins, nvars);
    for ii = 1:nbins
        arr2filter( (ii-1)*binsize+1:(ii-1)*binsize+burstsize(ii), :) = ...
            data_array(endpoint(ii)+1:endpoint(ii)+burstsize(ii), :);
    end
    arr2filter = reshape(arr2filter, binsize, nbins, nvars);

    %.. trap out the following situation which would result in no valid data at all.
    if nptsToSkip >= binsize
        error('blankingTimeSec is too long wrt burst time in median filtering operation.')
    end
    %.. if the number of points in a bin happens to be less than the number of blanking 
    %.. points to be deleted then the associated time value will be nan. Replace these
    %.. values, if they occur, with median filtered time values from the unblanked data.
    tf_time = strcmpi(var1D_filtered_new_names, 'time');  % will only be one time variable
    unblanked_median_time = nanmedian(arr2filter(:, :, tf_time), 1);
    
    arr2filter(1:nptsToSkip, :, :) = [];  % exclude data to be blanked at burst starts
    arr2filter = nanmedian(arr2filter, 1);
    filteredArray(1:nbins, 1:nvars) = arr2filter;
    %.. find locations of nans in the median filtered blanked time record
    tf_timeNans = isnan(filteredArray(:, tf_time));
    %.. replace them with values from the median filtered unblanked time record
    filteredArray(tf_timeNans, tf_time) = unblanked_median_time(tf_timeNans);
    
    %.. write out burst median filtered 1D variables
    for ii = 1:nvars
        ncwrite(outfilename, var1D_filtered_new_names{ii}, filteredArray(:, ii));
    end
    clearvars arr2filter filteredArray

    %.. MEDIAN FILTER 2D VARIABLES
    %.. process the 2D data variables one at a time
    nvars = length(idx_var2D);
    for jj = 1:nvars
        arr2filter = NaN(binsize*nbins, nwvl_acs);
        varname = schema.Variables(idx_var2D(jj)).Name;
        %.. transpose data array so that time changes by row number
        data_array = ncread(infilename, varname)';
        %.. change infinities to nan
        data_array(isinf(data_array(:))) = nan;
        for ii = 1:nbins
            arr2filter( (ii-1)*binsize+1:(ii-1)*binsize+burstsize(ii), :) = ...
                data_array(endpoint(ii)+1:endpoint(ii)+burstsize(ii), :);
        end
        arr2filter = reshape(arr2filter, binsize, nbins, nwvl_acs);
        arr2filter(1:nptsToSkip, :, :) = [];  % exclude data at start of the burst
        arr2filter = nanmedian(arr2filter, 1);
        filteredArray(1:nbins, 1:nwvl_acs) = arr2filter;
        %.. write out burst median filtered 2D variable
        ncwrite(outfilename, varname, filteredArray');
    end
end  % tf_burstMedian

%.. this section is independent of whether burstMedian filtering was applied.
%****************************************************************************
if tf_addSpectral1D
    %************************************************************************
    %.. CREATE FINAL SPECTRAL 1D DATA PRODUCTS BY INTERPOLATION 
    %************************************************************************
    %.. .. infinities have been replaced by nans in outfile
    %.. .. wavelength changes by row
    absorption = ncread(outfilename, 'optical_absorption');  % 2D array
    beam_c     = ncread(outfilename,  'beam_attenuation');   % 2D array
    %.. read in wavelengths
    wvl_a = ncread(outfilename, 'wavelength_a');
    wvl_c = ncread(outfilename, 'wavelength_c');
    %.. some netcdf files downloaded from OOINET can have wavelength values of
    %.. Nan added as padding to the end of the data vector. eliminate these if
    %.. they exist so that the interpolation to calculate the spectral 1D
    %.. products can be executed.
    tf_isnan = isnan(wvl_a);
    wvl_a(tf_isnan) = [];
    wvl_c(tf_isnan) = [];
    absorption(tf_isnan, :) = [];
    beam_c(tf_isnan, :)     = [];
    %.. interpolate to get values at the 1D wavelengths
    abs_1D    = interp1(wvl_a, absorption, wavelengths_1D');
    beam_c_1D = interp1(wvl_c, beam_c,     wavelengths_1D');
    
    %.. the 2D arrays are read in as [wavelength x time];
    %.. the 1D variables will be written out as column vectors (may not matter)
    for ii = 1:nwvl_1D
        ncwrite(outfilename, abs_1D_varname{ii},     abs_1D(ii, :)' );
        ncwrite(outfilename, beam_c_1D_varname{ii},  beam_c_1D(ii, :)');
    end
    %************************************************************************
    %.. END:  CREATE FINAL SPECTRAL 1D DATA PRODUCTS BY INTERPOLATION
    %************************************************************************
end  % tf_addSpectral1D
%****************************************************************************
