#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, load_gc_thredds, get_deployment_dates, \
    get_vocabulary, update_dataset, ENCODINGS
from ooi_data_explorations.profilers import create_profile_id
from ooi_data_explorations.qartod.qc_processing import parse_qc


def ctdpf_wfp(ds):
    """
    Takes CTD data recorded by the WFP and lightly re-works it into a xarray
    dataset for further work by the user. Primary task is dropping a
    variable, renaming three raw variables to be explicit and converting the
    original OOI QC variables from bitmaps to a QARTOD style flag.
    Additionally, adds a variable for the profile number.

    :param ds: initial CTDPF data set downloaded from OOI
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == time, redundant so can remove
    ds = ds.reset_coords()
    drop_vars = ['internal_timestamp']
    for v in drop_vars:
        if v in ds.variables:
            ds = ds.drop(v)

    # rename some of the variables for better clarity
    rename = {
        'conductivity': 'raw_sea_water_electrical_conductivity',
        'temperature': 'raw_sea_water_temperature',
        'pressure': 'raw_sea_water_pressure',
    }
    for key, value in rename.items():
        if key in ds.variables:
            ds = ds.rename({key: value})
            ds[value].attrs['ooinet_variable_name'] = key

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    # add a variable for the profile number
    ds = create_profile_id(ds)

    return ds


def ctdpf_cspp(ds):
    """
    Takes CTD data recorded by the CSPP controller and lightly re-works it into
    a xarray dataset for further work by the user. Primary task is dropping two
    redundant variables and converting the original OOI QC variables from
    bitmaps to a QARTOD style flag. Additionally, adds a variable for the
    profile number.

    :param ds: initial CTDPF data set downloaded from OOI
    :return: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == time, redundant so can remove
    #   profiler_timestamp == internal_timestamp, redundant so can remove
    drop_vars = ['internal_timestamp', 'profiler_timestamp']
    for v in drop_vars:
        if v in ds.variables:
            ds = ds.drop(v)

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    # add a variable for the profile number
    ds = create_profile_id(ds)

    return ds


def main(argv=None):
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor
    method = args.method
    stream = args.stream
    deploy = args.deploy
    start = args.start
    stop = args.stop

    # check if we are specifying a deployment or a specific date and time range
    if not deploy or (start and stop):
        return SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')

    # if we are specifying a deployment number, then get the data from the Gold Copy THREDDS server
    if deploy:
        # download the data for the deployment
        ctdpf = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*CTDPF.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not ctdpf:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, deployment %d.' % (site, node, sensor, method,
                                                                                    stream, deploy))
            raise SystemExit(exit_text)
    else:
        # otherwise, request the data for download from OOINet via the M2M API using the specified dates
        r = m2m_request(site, node, sensor, method, stream, start, stop)
        if not r:
            exit_text = ('Request failed for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                  stream, start, stop))
            raise SystemExit(exit_text)

        # Valid M2M request, start downloading the data
        ctdpf = m2m_collect(r, '.*CTDPF.*\\.nc$')

        # check to see if we downloaded any data
        if not ctdpf:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize the data
    if node == 'SP001':
        # this CTDPF is part of a CSPP
        ctdpf = ctdpf_cspp(ctdpf)
    else:
        # this CTDPF is on a WFP
        ctdpf = ctdpf_wfp(ctdpf)

    vocab = get_vocabulary(site, node, sensor)[0]
    ctdpf = update_dataset(ctdpf, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    ctdpf.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
