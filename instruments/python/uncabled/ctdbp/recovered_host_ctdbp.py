#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import yaml

from instruments.python.common import list_deployments, deployment_dates, get_vocabulary, m2m_request, m2m_collect, \
    update_dataset
from instruments.python.uncabled.ctdbp.request_ctdbp import ctdbp_datalogger

CONFIG = yaml.safe_load(open('instruments\\python\\config.yaml'))


def main():
    # Setup needed parameters for the request, the user would need to vary these to suit their own needs and
    # sites/instruments of interest. Site, node, sensor, stream and delivery method names can be obtained from the
    # Ocean Observatories Initiative web site. The last two will set path and naming conventions to save the data
    # to the local disk
    site = 'CE02SHSM'           # OOI Net site designator
    node = 'RID27'              # OOI Net node designator
    sensor = '03-CTDBPC000'     # OOI Net sensor designator
    stream = 'ctdbp_cdef_dcl_instrument_recovered'  # OOI Net stream name
    method = 'recovered_host'   # OOI Net data delivery method
    level = 'nsif'              # local directory name, level below site
    instrmt = 'ctdbp'           # local directory name, instrument below level

    # We are after recovered host data. Determine list of deployments and use the first deployment to determine
    # the start and end dates for our request.
    vocab = get_vocabulary(site, node, sensor)[0]
    deployments = list_deployments(site, node, sensor)
    deploy = deployments[0]
    start, stop = deployment_dates(site, node, sensor, deploy)

    # request and download the data
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    ctdbp = m2m_collect(r, '.*ctdbp.*\\.nc$')
    ctdbp = ctdbp.where(ctdbp.deployment == deploy, drop=True)  # limit to the deployment of interest

    # clean-up and reorganize
    ctdbp = ctdbp_datalogger(ctdbp, burst=True)
    ctdbp = update_dataset(ctdbp, vocab['maxdepth'])

    # save the data
    out_path = CONFIG['base_dir']['m2m_base'] + ('\\%s\\%s\\%s' % (site.lower(), level, instrmt))
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    out_file = ('%s.%s.%s.deploy%02d.%s.%s.nc' % (site.lower(), level, instrmt, deploy, method, stream))
    nc_out = os.path.join(out_path, out_file)

    ctdbp.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='netcdf4')


if __name__ == '__main__':
    main()