function revamp_OOINet_spkir_netcdf_files(infile, outfile)
%.. 08-dec-2018 desiderio 
%.. .. adapted from make_uframe_spkir_netcdf_test_files_1D_and_2D.m:
%.. .. (1) median filtering is removed
%.. .. (2) only one file, with both the 1D and 2D data products, is written out
%
%.. infile  = uframe SPKIR netcdf filename
%.. outfile = name of output netcdf file
%.. 20-dec-2019 desiderio. renamed function to revamp_OOINet_spkir_netcdf_files.m
%.. 19-jun-2020 desiderio. added 'pressure' to the vars2keep array. 

%.. operates on SPKIR (Satlantic OCR-507) UFrame netcdf files
%.. (1) deletes extraneous dimensions
%.. (2) deletes extraneous variables
%.. (3) replaces 'obs' dimension with 'time' to make time a coordinate variable
%.. (4) changes name of 'spectra' dimension to 'wavelength'
%.. (5) renames 'channel_array' variable to 'channel_counts'
%.. (6) renames 'spkir_abj_cspp_downwelling_vector' -> 'irradiance'
%.. (7) adds wavelength variable using the nominal OCR-507 wavelength values
%.. (8) renames units of irradiance variable to 'uW cm-2 nm-1'
%.. (9) RETAINS coordinate attribute of irradiance variable (previously deleted)
%.. (A) change 'irradiance' ancillary_variables attribute from
%..     'channel_array' -> 'channel_counts'
%.. (B) separates out the wavelength data into 7 channels (counts, irradiance)
%.. (C) writes out 1 netcdf file with both 1D and 2D data variables

%.. NOTE
%.. Matlab 2018b does not support the NC_STRING netcdf datatype. If this is
%.. encountered this code will throw a runtime execution error. To date the
%.. uframe netcdf files encountered have not used this datatype and instead
%.. have used NC_CHAR.

%************************************************************************
%***** CONDITION DIMENSIONS AND VARIABLES FROM INPUT NETCDF FILE ********
%************************************************************************
nominal_ocr507_wavelengths = [412 443 490 510 555 620 683];
nwvl = length(nominal_ocr507_wavelengths);

%.. old dimension names to keep from infile
dims2keep = {'obs' 'spectra'};
%.. these will be renamed as:
old_abscissa_dim_name = 'obs';
new_abscissa_dim_name = 'time';
%
old_wavelength_dim_name = 'spectra';
new_wavelength_dim_name = 'wavelength';

%.. variables (uframe names) to keep from infile. note that while 'spectra'
%.. is specified as a dimension, it is not a variable.
%.. .. for generality, three variable names are listed which have not been
%.. .. in nsif files but have been present at various times in cspp files:
%.. .. they are 'int_ctd_pressure', 'pressure_depth', and 'pressure'.
%.. 'spkir_abj_cspp_downwelling_vector' is used in both cspp and nsif files.
vars2keep = {
    'channel_array'
    'int_ctd_pressure'
    'internal_temperature'
    'lat'
    'lon'
    'pressure'
    'pressure_depth'
    'spkir_abj_cspp_downwelling_vector'
    'time'
    };
%.. data products to be renamed in outfile
oldRawProductName = 'channel_array';
rawProductName    = 'channel_counts';  % units of counts
%
oldDataProductName = 'spkir_abj_cspp_downwelling_vector';
dataProductName    = 'irradiance';     % units of uW cm-2 nm-1

schema = ncinfo(infile);

%.. delete extraneous dimensions
alldims = {schema.Dimensions.Name};
schema.Dimensions(~ismember(alldims, dims2keep)) = [];

%.. delete the variables not listed in vars2keep
allvars = {schema.Variables.Name};
schema.Variables(~ismember(allvars, vars2keep)) = [];

%.. nsif and cspp have different variable lists. to handle both cases when
%.. writing to outfile, get actual listing before schema is modified.
old_varname = {schema.Variables.Name};

%.. change the abscissa and wavelength dimension names from old to new
tf = strcmp({schema.Dimensions.Name}, old_abscissa_dim_name);
schema.Dimensions(tf).Name = new_abscissa_dim_name;
tf = strcmp({schema.Dimensions.Name}, old_wavelength_dim_name);
schema.Dimensions(tf).Name = new_wavelength_dim_name;

%.. change the variables' dimension assignations from old to new
for ii = 1:length(schema.Variables)
    for jj = 1:length(schema.Variables(ii).Dimensions)
        if strcmp(schema.Variables(ii).Dimensions(jj).Name, ...
                old_abscissa_dim_name)
            schema.Variables(ii).Dimensions(jj).Name = new_abscissa_dim_name;
        elseif strcmp(schema.Variables(ii).Dimensions(jj).Name, ...
                old_wavelength_dim_name)
            schema.Variables(ii).Dimensions(jj).Name = new_wavelength_dim_name;
        end
    end
end

%.. change the data product names 
tf = strcmp({schema.Variables.Name}, oldRawProductName);
schema.Variables(tf).Name = rawProductName;
tf = strcmp({schema.Variables.Name}, oldDataProductName);
schema.Variables(tf).Name = dataProductName;
%.. and make a new list documenting the changed names for output
new_varname = strrep(old_varname, oldRawProductName, rawProductName);
new_varname = strrep(new_varname, oldDataProductName, dataProductName);

%.. add another structure to schema containing wavelength-as-a-variable info;
%.. .. values taken from OMS-ERDDAP SPKIR netcdf file.
%.. construct inner structures first:
%.. .. Dimensions
dim.Name      = 'wavelength';
dim.Length    = nwvl;
dim.Unlimited = 0;
%.. .. Attributes
att = struct('Name' , {'units' 'standard_name'        'long_name'}, ...
             'Value', {'nm'    'radiation_wavelength' 'Radiation Wavelength'});
%.. structure for the wavelength variable
wvl.Name         = 'wavelength';
wvl.Dimensions   = dim;
wvl.Size         = nwvl;
wvl.Datatype     = 'int32';
wvl.Attributes   = att;
wvl.ChunkSize    = [];
wvl.FillValue    = -2147483647;
wvl.DeflateLevel = [];
wvl.Shuffle      = false;
%.. write as another Variable array element into schema
schema.Variables(end+1) = wvl;

%.. find units for the irradiance data product which has been read in to
%.. matlab as 'UNSUPPORTED DATATYPE' because Uframe used the greek letter
%.. {micro}, which matlab does not support, instead of a lower case 'u'. 
tf_name = strcmp({schema.Variables.Name}, dataProductName);
tf_units = strcmp({schema.Variables(tf_name).Attributes.Name}, 'units');
schema.Variables(tf_name).Attributes(tf_units).Value = 'uW cm-2 nm-1';
%.. also change the name of the raw data product name in:
tf_ancillvar = strcmp({schema.Variables(tf_name).Attributes.Name}, ...
    'ancillary_variables');
schema.Variables(tf_name).Attributes(tf_ancillvar).Value = rawProductName;

%.. this version: keep
% %.. DELETE the (auxiliary) coordinate attribute to this variable,
% %.. in case it may be slowing things up; note that CWingard's OMS-ERDDAP
% %.. 'auto'-generated spkir netcdf files do *not* have a coordinate
% %.. attribute assigned to the downwelling irradiance data product.
% tf_coordinates = strcmp({schema.Variables(tf_name).Attributes.Name}, ...
%     'coordinates');
% schema.Variables(tf_name).Attributes(tf_coordinates) = [];

%************************************************************************
%.. ADD NETCDF INFO FOR 1D COUNTS AND Ed VARIABLES TO SCHEMA
%************************************************************************

%.. construct spectral 1D variable names
cts_1D_varname(1:nwvl) = {''};
Ed_1D_varname(1:nwvl)  = {''};
for ii = 1:nwvl
    wave_char = [num2str(nominal_ocr507_wavelengths(ii)) 'nm'];
    cts_1D_varname{ii} = ['counts_' wave_char];
    Ed_1D_varname{ii}  = ['Ed_' wave_char];
end

%.. add netcdf schema specifications for nwvl raw counts 1D variables
%.. .. isolate 2D channel_counts Variables structure to serve as a template
tf = strcmp({schema.Variables.Name}, rawProductName);
cts_template = schema.Variables(tf);
%.. reset dimensions from 2D to 1D by eliminating wavelength dimension
tf = strcmp({cts_template.Dimensions.Name}, 'wavelength');
cts_template.Dimensions(tf) = [];
%.. reset Size 2D -> 1D
cts_template.Size(cts_template.Size==nwvl) = [];
%.. reset ChunkSize 2D -> 1D
cts_template.ChunkSize(cts_template.ChunkSize==nwvl) = [];
%.. add the nwvl variables to the schema structure
for ii = 1:nwvl
    cts_template.Name = cts_1D_varname{ii};
    schema.Variables(end+1) = cts_template;
end

%.. add netcdf schema specifications for nwvl irradiance 1D variables
%.. .. isolate 2D irradiance Variables structure to serve as a template
tf = strcmp({schema.Variables.Name}, dataProductName);
irr_template = schema.Variables(tf);
%.. reset dimensions from 2D to 1D by eliminating wavelength dimension
tf = strcmp({irr_template.Dimensions.Name}, 'wavelength');
irr_template.Dimensions(tf) = [];
%.. reset Size 2D -> 1D
irr_template.Size(irr_template.Size==nwvl) = [];
%.. reset ChunkSize 2D -> 1D
irr_template.ChunkSize(irr_template.ChunkSize==nwvl) = [];
%.. to reassign ancillary variables attribute to corresponding counts channel:
tf_ancillary = strcmp({irr_template.Attributes.Value}, rawProductName);
%.. add the nwvl variables to the schema structure
for ii = 1:nwvl
    irr_template.Name = Ed_1D_varname{ii};
    irr_template.Attributes(tf_ancillary).Value = cts_1D_varname{ii};
    schema.Variables(end+1) = irr_template;
end


%.. WRITE TO OUTFILE
if ~strcmp(outfile(end-2:end), '.nc')
    outfile = [outfile '.nc'];
end
%.. MAKE SURE THAT THIS FILE DOES NOT ALREADY EXIST!
%.. by design the ncwriteschema function *appends* schema to the file,
%.. if it exists; in this code we want a completely new file.
if isfile(outfile)
    error('New file exists; it MUST be removed for code to work properly.');
end
%.. create empty netcdf file with new schema.
ncwriteschema(outfile, schema);

%..write out old data first
for ii = 1:length(new_varname)
    tmp = ncread(infile, old_varname{ii});
    ncwrite(outfile, new_varname{ii}, tmp);
end
%.. and write in the SPKIR wavelength values to the wavelength variable
ncwrite(outfile, 'wavelength', nominal_ocr507_wavelengths);

%.. now read in 2D data so that spectral 1D data can be written
counts = ncread(infile, oldRawProductName);       % 2D array
Ed     = ncread(infile, oldDataProductName);      % 2D array
%.. the 2D arrays are read in as [wavelength x time];
%.. the 1D variables will be written out as column vectors (may not matter)
for ii = 1:nwvl
    ncwrite(outfile, cts_1D_varname{ii}, counts(ii, :)' );
    ncwrite(outfile, Ed_1D_varname{ii},  Ed(ii, :)'     );
end
