#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import m2m_request, m2m_collect, load_gc_thredds, get_deployment_dates, \
    update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_vel3d import vel3d_datalogger

# Setup needed parameters for the request, the user would need to vary
# these to suit their own needs and sites/instruments of interest. Site,
# node, sensor, stream and delivery method names can be obtained from the
# Ocean Observatories Initiative web site. The last two parameters (level
# and instrmt) will set path and naming conventions to save the data to the
# local disk.
site = 'CE01ISSM'           # OOI Net site designator
node = 'MFD35'              # OOI Net node designator
sensor = '01-VEL3DD000'     # OOI Net sensor designator
method = 'telemetered'      # OOI Net data delivery method
level = 'mfn'              # local directory name, level below site
instrmt = 'vel3d'           # local directory name, instrument below level

# download the data for deployment 17 from the Gold Copy THREDDS catalog
tag = 'deployment0017.*VEL3D.*\\.nc$'
# -- the velocity data is available in the GC THREDDS catalog
velocity = load_gc_thredds(site, node, sensor, method, 'vel3d_cd_dcl_velocity_data', tag)
# -- the additional system and header data are available in the M2M system because they are not
# -- considered "science" data, but rather "engineering" data which is not included in the GC THREDDS catalog.
start, stop = get_deployment_dates(site, node, sensor, 17)
m = m2m_request(site, node, sensor, method, 'vel3d_cd_dcl_system_data', start, stop)
system = m2m_collect(m, tag)  # provides heading, pitch, roll, and temperature data
m = m2m_request(site, node, sensor, method, 'vel3d_cd_dcl_data_header', start, stop)
header = m2m_collect(m, tag)  # provides noise amplitude measurements used to QC the data

# clean-up and reorganize
vel3d = vel3d_datalogger(velocity, system, header, burst=True)


# save the data
out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), level, instrmt)
out_path = os.path.abspath(out_path)
if not os.path.exists(out_path):
    os.makedirs(out_path)

out_file = ('%s.%s.%s.deploy17.%s.%s.nc' % (site.lower(), level, instrmt, method, stream))
nc_out = os.path.join(out_path, out_file)

velpt.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)
