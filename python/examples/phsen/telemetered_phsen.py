#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import list_deployments, get_deployment_dates, get_vocabulary, m2m_request, \
    m2m_collect, update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_phsen import phsen_datalogger


def main():
    # Setup needed parameters for the request, the user would need to vary these to suit their own needs and
    # sites/instruments of interest. Site, node, sensor, stream and delivery method names can be obtained from the
    # Ocean Observatories Initiative web site. The last two will set path and naming conventions to save the data
    # to the local disk
    site = 'CE02SHSM'           # OOI Net site designator
    node = 'RID26'              # OOI Net node designator
    sensor = '06-PHSEND000'     # OOI Net sensor designator
    stream = 'phsen_abcdef_dcl_instrument'  # OOI Net stream name
    method = 'telemetered'      # OOI Net data delivery method
    level = 'nsif'              # local directory name, level below site
    instrmt = 'phsen'           # local directory name, instrument below level

    # We are after telemetered data. Determine list of deployments and use the last, presumably currently active,
    # deployment to determine the start and end dates for our request.
    vocab = get_vocabulary(site, node, sensor)[0]
    deployments = list_deployments(site, node, sensor)
    deploy = deployments[-1]
    start, stop = get_deployment_dates(site, node, sensor, deploy)

    # request and download the data
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    phsen = m2m_collect(r, '.*PHSEN.*\\.nc$')
    phsen = phsen.where(phsen.deployment == deploy, drop=True)  # limit to the deployment of interest

    # clean-up and reorganize
    phsen = phsen_datalogger(phsen)
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
