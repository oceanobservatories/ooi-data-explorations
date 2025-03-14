#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import get_vocabulary, load_gc_thredds, update_dataset, \
    CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_wavss import wavss_directional

# Setup needed parameters for the request, the user would need to vary these to suit their own needs and
# sites/instruments of interest. Site, node, sensor, stream and delivery method names can be obtained from the
# Ocean Observatories Initiative website. The last two parameters (level and instrmt) will set path and naming
# conventions to save the data to the local disk.
site = 'CE02SHSM'           # OOI Net site designator
node = 'SBD12'              # OOI Net node designator
sensor = '05-WAVSSA000'     # OOI Net sensor designator
stream = 'wavss_a_dcl_mean_directional_recovered'  # OOI Net stream name
method = 'recovered_host'   # OOI Net data delivery method
level = 'buoy'              # local directory name, level below site
instrmt = 'wavss'           # local directory name, instrument below level

# We are after recovered_host data from deployment 13, download the data from the OOI Gold Copy THREDDS catalog
vocab = get_vocabulary(site, node, sensor)[0]
deploy = 13
wavss = load_gc_thredds(site, node, sensor, method, stream, '.*deployment%04d.*WAVSS.*mean_directional.*\\.nc$' % deploy)

# clean-up and reorganize the WAVSS data set
wavss = wavss_directional(wavss)
wavss = update_dataset(wavss, vocab['maxdepth'])

# save the data -- utilize groups for the wavss and water datasets
out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), level, instrmt)
out_path = os.path.abspath(out_path)
if not os.path.exists(out_path):
    os.makedirs(out_path)

out_file = ('%s.%s.%s.deploy%02d.%s.%s.nc' % (site.lower(), level, instrmt, deploy, method, stream))
nc_out = os.path.join(out_path, out_file)
wavss.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)
