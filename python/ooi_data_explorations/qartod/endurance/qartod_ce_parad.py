#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Load the PARAD data from the uncabled, Coastal Endurance Profiler
    Mooring and Coastal Surface-Piercing Profilers (CSPP) and process the data
    to generate QARTOD Gross Range and Climatology test limits
"""
import dask
import numpy as np
import os
import pandas as pd
import pytz
import xarray as xr

from ftplib import FTP
from pysolar.solar import get_altitude
from pysolar.radiation import get_radiation_direct

from ooi_data_explorations.common import get_annotations, get_vocabulary, list_deployments, get_deployment_dates, \
    get_sensor_information
from ooi_data_explorations.qartod.qc_processing import process_gross_range, process_climatology, woa_standard_bins, \
    inputs, ANNO_HEADER, CLM_HEADER, GR_HEADER


# noinspection PyTypeChecker
def generate_qartod(site, node, sensor):
    """
    blah, blah, blah.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :return gr_lookup: CSV formatted strings to save to a csv file for the
        QARTOD gross range lookup tables.
    :return clm_lookup: CSV formatted strings to save to a csv file for the
        QARTOD climatology lookup tables.
    :return clm_table: CSV formatted strings to save to a csv file for the
        QARTOD climatology range tables.
    """
    # get the current system annotations for the sensor
    annotations = get_annotations(site, node, sensor)
    annotations = pd.DataFrame(annotations)
    if not annotations.empty:
        annotations = annotations.drop(columns=['@class'])
        annotations['beginDate'] = pd.to_datetime(annotations.beginDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')
        annotations['endDate'] = pd.to_datetime(annotations.endDT, unit='ms').dt.strftime('%Y-%m-%dT%H:%M:%S')

    # get the min and max depth for the site-node-sensor from the M2M vocabulary
    vocab = get_vocabulary(site, node, sensor)[0]
    min_depth = vocab['mindepth']
    max_depth = vocab['maxdepth']

    # gather the latitude and longitude data for all the deployments and calculate the mean location
    deployments = list_deployments(site, node, sensor)
    latitude = []
    longitude = []
    for deploy in deployments:
        metadata = get_sensor_information(site, node, sensor, deploy)[0]
        latitude.append(metadata['location']['latitude'])
        longitude.append(metadata['location']['longitude'])

    latitude = np.mean(latitude)
    longitude = np.mean(longitude)

    # set up, if not already created a directory to save the Coast Watch data (saves reprocessing the data)
    kd_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/coastwatch/kdpar')
    kd_path = os.path.abspath(kd_path)
    if not os.path.exists(kd_path):
        os.makedirs(kd_path)

    kd_file = '%s.noaa.viirs.monthly.sq.kdpar.nc' % site.lower()
    kd_out = os.path.join(kd_path, kd_file)
    if os.path.isfile(kd_out):
        kd = xr.load_dataset(kd_out)
    else:
        # Data obtained from the "NOAA MSL12 Ocean Color - Science Quality - VIIRS SNPP" satellite data products
        # page, downloading the L3 monthly KdPAR data from the FTP server (ERDDAP and THREDDS servers proved too
        # unstable to rely on). See the data products web page for more information:
        #   https://coastwatch.noaa.gov/cw/satellite-data-products/ocean-color/science-quality/viirs-snpp.html
        # check to see if the data has been downloaded, if not do so first
        if not os.path.exists(os.path.join(kd_path, 'raw')):
            os.makedirs(os.path.join(kd_path, 'raw'))
            ftp_server = 'ftp.star.nesdis.noaa.gov'
            ftp_server_path = '/pub/socd1/mecb/coastwatch/viirs/science/L3/global/kd/monthly/WW00/'
            ftp = FTP(ftp_server)
            ftp.login(user='anonymous')
            ftp.cwd(ftp_server_path)
            files = ftp.nlst('*kdpar.nc')
            for file in files:
                print('downloading %s' % os.path.join(kd_path, 'raw', os.path.basename(file)))
                with open(os.path.join(kd_path, 'raw', os.path.basename(file)), 'wb') as f:
                    ftp.retrbinary('RETR %s' % file, f.write)

            ftp.quit()

        # load the Kd(PAR) data into an xarray dataset
        with dask.config.set(**{'array.slicing.split_large_chunks': False}):
            kd = xr.open_mfdataset(os.path.join(kd_path, 'raw/*.nc'), combine='nested', concat_dim='time',
                                   engine='netcdf4')

        # clean up the data set and use a simple bounding box to limit the extent to the site of interest
        kd = kd.drop_vars(['coord_ref', 'palette'])
        kd = kd.where((kd['lat'] >= latitude - 0.09375) & (kd['lat'] <= latitude + 0.09375), drop=True)
        kd = kd.where((kd['lon'] >= longitude - 0.09375) & (kd['lon'] <= longitude + 0.09375), drop=True)

        kd = kd.mean(dim=['altitude', 'lat', 'lon'], keep_attrs=True)
        kd = kd.sortby('time')

        kd['kd_par'] = kd.kd_par.compute()  # convert dask array to standard
        kd.to_netcdf(kd_out, mode='w', format='NETCDF4', engine='h5netcdf')

    # calculate clear sky irradiance at solar noon for this site using the Kd(PAR) time record
    date = pd.to_datetime(kd.time).to_pydatetime()
    surface = []
    for d in date:
        dt = d.replace(tzinfo=pytz.timezone('US/Pacific'))
        altitude = get_altitude(latitude, longitude, dt)
        # PAR is approximately 50% of the shortwave radiation, and we need to convert from W/m^2 to umol/m^2/s
        surface.append(get_radiation_direct(dt, altitude) * 0.5 / 0.21739130434)

    # create a 2D array with Ed(PAR) estimated as a function of depth from the satellite Kd(PAR) values and
    # model estimates of clear-sky irradiance
    depths = np.arange(min_depth, max_depth + 0.125, 0.125)
    ed = np.zeros([len(date), len(depths)])
    for i in range(len(date)):
        ed[i, :] = surface[i] * np.exp(-kd.kd_par.values[i] * depths)

    # convert to an xarray dataset
    ed = xr.Dataset({
        'Ed': (['time', 'depth'], ed),
    }, coords={'time': kd.time.values, 'depth': depths})

    # set the parameters and the gross range limits
    parameters = ['Ed']
    limits = [0, 5000]

    # create the initial gross range entry
    quantile = ed['depth'].quantile(0.01).values     # upper 1% of the depth array
    sub = ed.where(ed.depth <= quantile, drop=True)  # limit gross range estimation to near-surface values
    sub = sub.max(dim=['depth'], keep_attrs=True)
    gr_lookup = process_gross_range(sub, parameters, limits, site=site,
                                    node=node, sensor=sensor, stream='parad_replace_me')

    # add the stream name and the source comment
    gr_lookup['notes'] = ('User range modeled from data collected by the NOAA VIIRS satellite and estimates of '
                          'clear sky irradiance from the pysolar package.')

    # create the depth bins for a depth-based climatology
    depth_bins = woa_standard_bins()
    m = depth_bins[:, 1] <= max_depth
    depth_bins = depth_bins[m, :]
    if site == 'CE09OSPM':
        depth_bins = depth_bins[6:, :]  # WFP doesn't go above 30 m

    # create and format the climatology lookups and tables for the data
    clm_lookup, clm_table = process_climatology(ed, parameters, limits, depth_bins=depth_bins,
                                                site=site, node=node, sensor=sensor,
                                                stream='parad_replace_me', fixed_lower=True)

    return annotations, gr_lookup, clm_lookup, clm_table


def main(argv=None):
    """
    Download the PARAD data from the Gold Copy THREDDS server and create the
    QARTOD gross range and climatology test lookup tables.
    """
    # set up the input arguments
    args = inputs(argv)
    site = args.site
    node = args.node
    sensor = args.sensor

    # create the QARTOD gross range and climatology lookup values and tables
    annotations, gr_lookup, clm_lookup, clm_table = generate_qartod(site, node, sensor)

    # save the downloaded annotations and qartod lookups and tables
    out_path = os.path.join(os.path.expanduser('~'), 'ooidata/qartod/parad')
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # save the annotations to a csv file for further processing
    anno_csv = '-'.join([site, node, sensor]) + '.quality_annotations.csv'
    annotations.to_csv(os.path.join(out_path, anno_csv), index=False, columns=ANNO_HEADER)

    # save the gross range values to a csv for further processing
    gr_csv = '-'.join([site, node, sensor]) + '.gross_range.csv'
    gr_lookup.to_csv(os.path.join(out_path, gr_csv), index=False, columns=GR_HEADER)

    # save the climatology values and table to a csv for further processing
    clm_csv = '-'.join([site, node, sensor]) + '.climatology.csv'
    clm_lookup.to_csv(os.path.join(out_path, clm_csv), index=False, columns=CLM_HEADER)
    parameters = ['PAR']
    for i in range(len(parameters)):
        tbl = '-'.join([site, node, sensor, parameters[i]]) + '.csv'
        with open(os.path.join(out_path, tbl), 'w') as clm:
            clm.write(clm_table[i])


if __name__ == '__main__':
    main()
