function variables = M2M_OPTAA(netcdfFilenames, tf_addDiscreteWavelengthTimeSeries)

netcdfFilenames = string(netcdfFilenames);  % so that character vector will work
ooifile = netcdfFilenames(contains(netcdfFilenames, 'OPTAA'));
gridded_file = strcat('processed_', ooifile);
for i = 1:length(ooifile)
    revamp_OOINet_optaa_netcdf_files(ooifile{i},'dummy.nc',[],tf_addDiscreteWavelengthTimeSeries);
    grid_optaa_netcdf_wavelength('dummy.nc',gridded_file{i});
    delete('dummy.nc')
end

%Read in wavelengths
variables.wavelength = ncread(gridded_file{1},'wavelength');

%Cat 2D attenuation and 2D absorption data from different files
mtime=[];attenuation_spectra_2d=[];absorption_spectra_2d=[];
for i = 1:length(ooifile)
    time=ncread(gridded_file{i},'time');
    time=datenum(1900,1,1,0,0,0)+time/60/60/24;
    mtime(length(mtime)+1:length(mtime)+length(time)) = time;
    
    beam_attenuation=ncread(gridded_file{i},'beam_attenuation');
    optical_absorption=ncread(gridded_file{i},'optical_absorption');
    
    %.. time increases as column number increases
    attenuation_spectra_2d = horzcat(attenuation_spectra_2d, beam_attenuation);
    absorption_spectra_2d = horzcat(absorption_spectra_2d, optical_absorption);
end
variables.mtime = mtime;
variables.attenuation_spectra_2d = attenuation_spectra_2d;
variables.absorption_spectra_2d = absorption_spectra_2d;

%Cat 1D spectral time series data
if tf_addDiscreteWavelengthTimeSeries
    nc_struct = ncinfo(gridded_file{1});
    vars = {nc_struct.Variables.Name}';
    vars = vars(contains(vars, {'abs_' 'beam_c_'}));
    for ii = 1:length(vars)
        varII=[];
        for jj = 1:length(ooifile)
            tmpJJ=ncread(gridded_file{jj}, vars{ii});
            varII(length(varII)+1:length(varII)+length(tmpJJ)) = tmpJJ;
        end
        variables.(vars{ii}) = varII;
    end
end    
    
end