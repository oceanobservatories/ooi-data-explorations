#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import list_deployments, get_vocabulary, load_gc_thredds, \
    update_dataset, CONFIG, ENCODINGS
from ooi_data_explorations.uncabled.process_nutnr import suna_cspp


def main():
    # Setup needed parameters for the request, the user would need to vary these to suit their own needs and
    # sites/instruments of interest. Site, node, sensor, stream and delivery method names can be obtained from the
    # Ocean Observatories Initiative website. The last parameter (instrmt) will set path and naming
    # conventions to save the data to the local disk.
    site = 'CE02SHSP'           # OOI Net site designator
    node = 'SP001'              # OOI Net node designator
    sensor = '05-NUTNRJ000'     # OOI Net sensor designator
    stream = 'nutnr_j_cspp_instrument_recovered'  # OOI Net stream name
    method = 'recovered_cspp'   # OOI Net data delivery method
    instrmt = 'nutnr'           # local directory name, instrument below site

    # We are after the recovered data. Determine list of deployments and use data from one of the earlier deployments
    vocab = get_vocabulary(site, node, sensor)[0]
    deployments = list_deployments(site, node, sensor)
    deploy = deployments[-4]

    # download the data from the Gold Copy THREDDS server
    nutnr = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*NUTNR.*\\.nc$' % deploy))

    # clean-up and reorganize
    nutnr = suna_cspp(nutnr)
    nutnr = update_dataset(nutnr, vocab['maxdepth'])

    # save the data
    out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), instrmt)
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    out_file = ('%s.%s.deploy%02d.%s.%s.nc' % (site.lower(), instrmt, deploy, method, stream))
    nc_out = os.path.join(out_path, out_file)

    nutnr.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
