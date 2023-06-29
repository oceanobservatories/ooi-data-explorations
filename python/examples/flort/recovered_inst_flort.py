#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import list_deployments, get_vocabulary, load_gc_thredds, \
    update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_flort import flort_instrument


def main():
    # Setup needed parameters for the request, the user would need to vary
    # these to suit their own needs and sites/instruments of interest. Site,
    # node, sensor, stream and delivery method names can be obtained from the
    # Ocean Observatories Initiative website. The last two parameters (level
    # and instrmt) will set path and naming conventions to save the data to the
    # local disk.
    site = 'CE01ISSM'           # OOI Net site designator
    node = 'SBD17'              # OOI Net node designator
    sensor = '06-FLORTD000'     # OOI Net sensor designator
    stream = 'flort_sample'     # OOI Net stream name
    method = 'recovered_inst'   # OOI Net data delivery method
    level = 'buoy'              # local directory name, level below site
    instrmt = 'flort'           # local directory name, instrument below level

    # We are after recovered instrument data, let's use the 6th deployment cause why not?
    vocab = get_vocabulary(site, node, sensor)[0]
    deployments = list_deployments(site, node, sensor)
    deploy = deployments[5]

    # download the data
    tag = 'deployment{:04g}.*FLORT.*\\.nc$'.format(deploy)
    flort = load_gc_thredds(site, node, sensor, method, stream, tag)

    # clean-up and reorganize
    flort = flort_instrument(flort)
    flort = update_dataset(flort, vocab['maxdepth'])

    # save the data
    out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), level, instrmt)
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    out_file = ('%s.%s.%s.deploy%02d.%s.%s.nc' % (site.lower(), level, instrmt, deploy, method, stream))
    nc_out = os.path.join(out_path, out_file)

    flort.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
