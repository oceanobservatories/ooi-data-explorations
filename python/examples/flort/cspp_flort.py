#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import list_deployments, get_vocabulary, load_gc_thredds, \
    update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_flort import flort_cspp


def main():
    # Setup needed parameters for the request, the user would need to vary these to suit their own needs and
    # sites/instruments of interest. Site, node, sensor, stream and delivery method names can be obtained from the
    # Ocean Observatories Initiative web site. The last two parameters (level and instrmt) will set path and naming
    # conventions to save the data to the local disk.
    site = 'CE02SHSP'           # OOI Net site designator
    node = 'SP001'              # OOI Net node designator
    sensor = '07-FLORTJ000'     # OOI Net sensor designator
    stream = 'flort_sample'     # OOI Net stream name
    method = 'recovered_cspp'   # OOI Net data delivery method
    instrmt = 'flort'           # local directory name, instrument below site

    # We are after the recovered data. Determine list of deployments and use data from one of the earlier deployments
    vocab = get_vocabulary(site, node, sensor)[0]
    deployments = list_deployments(site, node, sensor)
    deploy = deployments[-4]

    # download the data from the Gold Copy THREDDS server
    flort = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*FLORT.*\\.nc$' % deploy))

    # clean-up and reorganize
    flort = flort_cspp(flort)
    flort = update_dataset(flort, vocab['maxdepth'])

    # save the data
    out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), instrmt)
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    out_file = ('%s.%s.deploy%02d.%s.%s.nc' % (site.lower(), instrmt, deploy, method, stream))
    nc_out = os.path.join(out_path, out_file)

    flort.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
