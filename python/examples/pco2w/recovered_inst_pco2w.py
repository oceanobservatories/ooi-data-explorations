#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import list_deployments, get_deployment_dates, get_vocabulary, m2m_request, \
    m2m_collect, update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_pco2w import pco2w_instrument


def main():
    # Setup needed parameters for the request, the user would need to vary these to suit their own needs and
    # sites/instruments of interest. Site, node, sensor and stream names can be obtained from the Ocean Observatories
    # Initiative web site
    site = 'CE01ISSM'           # OOI Net site designator
    node = 'RID16'              # OOI Net node designator
    sensor = '05-PCO2WB000'     # OOI Net sensor designator
    stream = 'pco2w_abc_instrument'  # OOI Net stream name
    method = 'recovered_inst'   # OOI Net data delivery method
    level = 'nsif'              # local directory name, level below site
    instrmt = 'pco2w'           # local directory name, instrument below level

    # We are after recovered host data. Determine list of deployments and use a previous one to determine the
    # start and end dates for our request.
    vocab = get_vocabulary(site, node, sensor)[0]
    deployments = list_deployments(site, node, sensor)
    deploy = deployments[-3]
    start, stop = get_deployment_dates(site, node, sensor, deploy)

    # request and download the data
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    pco2w = m2m_collect(r, '^(?!.*blank).*PCO2W.*nc$')
    pco2w = pco2w.where(pco2w.deployment == deploy, drop=True)  # limit to the deployment of interest

    # clean-up and reorganize
    pco2w = pco2w_instrument(pco2w)
    pco2w = update_dataset(pco2w, vocab['maxdepth'])

    # save the data
    out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), level, instrmt)
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    out_file = ('%s.%s.%s.deploy%02d.%s.%s.nc' % (site.lower(), level, instrmt, deploy, method, stream))
    nc_out = os.path.join(out_path, out_file)

    pco2w.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
