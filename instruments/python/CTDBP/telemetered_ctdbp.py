#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import yaml

from instruments.python.common import list_deployments, deployment_dates, m2m_request, m2m_collect
from instruments.python.CTDBP.request_ctdbp import ctdbp_datalogger

CONFIG = yaml.safe_load(open('instruments\\python\\config.yaml'))


def main():
    # Setup needed parameters for the request, the user would need to vary these to suit their own needs and
    # sites/instruments of interest. Site, node, sensor and stream names can be obtained from the Ocean Observatories
    # Initiative web site
    site = 'CE02SHSM'
    node = 'RID27'
    sensor = '03-CTDBPC000'
    stream = 'ctdbp_cdef_dcl_instrument'

    # We are after telemetered data. In other words, we only want the data from the most recent, current (presumably)
    # deployment. Determine list of deployments and use the most recent to determine the start and end dates for our
    # request.
    deployments = list_deployments(site, node, sensor)
    start, stop = deployment_dates(site, node, sensor, deployments[-1])

    # request and download the data
    r = m2m_request(site, node, sensor, 'telemetered', stream, start, stop)
    ctdbp = m2m_collect(r, '.*CTDBP.*\\.nc')
    ctdbp = ctdbp.where(ctdbp.deployment == deployments[-1], drop=True)  # limit to the most recent deployment
    ctdbp = ctdbp_datalogger(ctdbp)

    # save the data
    out_path = CONFIG['base_dir']['m2m_base'] + '\\' + site.lower() + '\\nsif\\ctdbp'
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    out_file = site.lower() + '.nsif.ctdbp.deploy' + ('%02d.' % deployments[-1]) + 'telemetered.' + stream + '.nc'
    nc_out = os.path.abspath(os.path.join(out_path, out_file))

    ctdbp.to_netcdf(nc_out, mode='a', format='NETCDF4', engine='netcdf4')


if __name__ == '__main__':
    main()