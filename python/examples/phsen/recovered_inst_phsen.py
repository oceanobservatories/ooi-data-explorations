#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import list_deployments, get_deployment_dates, get_vocabulary, m2m_request, \
    m2m_collect, update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_phsen import phsen_instrument


def main():
    # Setup needed parameters for the request, the user would need to vary these to suit their own needs and
    # sites/instruments of interest. Site, node, sensor and stream names can be obtained from the Ocean Observatories
    # Initiative web site
    site = 'CE06ISSM'           # OOI Net site designator
    node = 'RID16'              # OOI Net node designator
    sensor = '06-PHSEND000'     # OOI Net sensor designator
    stream = 'phsen_abcdef_instrument'  # OOI Net stream name
    method = 'recovered_inst'   # OOI Net data delivery method
    level = 'nsif'              # local directory name, level below site
    instrmt = 'phsen'           # local directory name, instrument below level

    # We are after recovered instrument data. Determine list of deployments and use a recent one to determine
    # the start and end dates for our request.
    vocab = get_vocabulary(site, node, sensor)[0]
    deployments = list_deployments(site, node, sensor)
    deploy = deployments[-3]
    start, stop = get_deployment_dates(site, node, sensor, deploy)

    # request and download the data
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    phsen = m2m_collect(r, '.*phsen.*\\.nc$')
    phsen = phsen.where(phsen.deployment == deploy, drop=True)  # limit to the deployment of interest

    # clean-up and reorganize
    phsen = phsen_instrument(phsen)
    phsen = update_dataset(phsen, vocab['maxdepth'])

    # save the data
    out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), level, instrmt)
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    out_file = ('%s.%s.%s.deploy%02d.%s.%s.nc' % (site.lower(), level, instrmt, deploy, method, stream))
    nc_out = os.path.join(out_path, out_file)

    phsen.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
