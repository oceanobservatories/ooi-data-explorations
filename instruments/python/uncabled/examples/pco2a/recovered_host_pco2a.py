#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from instruments.python.common import list_deployments, get_deployment_dates, get_vocabulary, m2m_request, m2m_collect, \
    update_dataset, CONFIG, ENCODINGS
from instruments.python.uncabled.request_pco2a import pco2a_datalogger


def main():
    # Setup needed parameters for the request, the user would need to vary these to suit their own needs and
    # sites/instruments of interest. Site, node, sensor, stream and delivery method names can be obtained from the
    # Ocean Observatories Initiative web site. The last two parameters (level and instrmt) will set path and naming
    # conventions to save the data to the local disk.
    site = 'CE02SHSM'           # OOI Net site designator
    node = 'SBD12'              # OOI Net node designator
    sensor = '04-PCO2AA000'     # OOI Net sensor designator
    stream = 'pco2a_a_dcl_instrument_air_recovered'  # OOI Net stream name
    method = 'recovered_host'      # OOI Net data delivery method
    level = 'buoy'              # local directory name, level below site
    instrmt = 'pco2a'           # local directory name, instrument below level

    # We are after recovered_host data. Determine list of deployments and use the last, presumably currently active,
    # deployment to determine the start and end dates for our request.
    vocab = get_vocabulary(site, node, sensor)[0]
    deployments = list_deployments(site, node, sensor)
    deploy = deployments[-1]
    start, stop = get_deployment_dates(site, node, sensor, deploy)

    # request and download the data -- air measurements
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    air = m2m_collect(r, ('.*deployment%04d.*PCO2A.*air.*\\.nc$' % deploy))
    air = air.where(air.deployment == deploy, drop=True)  # limit to the deployment of interest

    # request and download the data -- water measurements
    r = m2m_request(site, node, sensor, method, 'pco2a_a_dcl_instrument_water_recovered', start, stop)
    water = m2m_collect(r, ('.*deployment%04d.*PCO2A.*water.*\\.nc$' % deploy))
    water = water.where(water.deployment == deploy, drop=True)  # limit to the deployment of interest

    # clean-up and reorganize the air and water datasets
    air = pco2a_datalogger(air, True)
    air = update_dataset(air, vocab['maxdepth'])
    water = pco2a_datalogger(water, True)
    water = update_dataset(water, vocab['maxdepth'])

    # save the data -- utilize groups for the air and water datasets
    out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), level, instrmt)
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    out_file = ('%s.%s.%s.deploy%02d.%s.%s.nc' % (site.lower(), level, instrmt, deploy, method, stream))
    nc_out = os.path.join(out_path, out_file)
    air.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS, group='air')
    water.to_netcdf(nc_out, mode='a', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS, group='water')


if __name__ == '__main__':
    main()
