#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os

from ooi_data_explorations.common import inputs, load_gc_thredds, m2m_collect, m2m_request, get_vocabulary, \
    update_dataset, ENCODINGS
from ooi_data_explorations.qartod.qc_processing import parse_qc

# load configuration settings
FILL_INT = -9999999
ATTRS = dict({
    # dataset attributes --> raw values
    'raw_irradiance_412': {
        'long_name': 'Raw Downwelling Irradiance 412 nm',
        'comment': ('Raw downwelling spectral irradiance at 412 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-412',
    },
    'raw_irradiance_444': {
        'long_name': 'Raw Downwelling Irradiance 444 nm',
        'comment': ('Raw downwelling spectral irradiance at 444 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-444',
    },
    'raw_irradiance_490': {
        'long_name': 'Raw Downwelling Irradiance 490 nm',
        'comment': ('Raw downwelling spectral irradiance at 490 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-490',
    },
    'raw_irradiance_510': {
        'long_name': 'Raw Downwelling Irradiance 510 nm',
        'comment': ('Raw downwelling spectral irradiance at 510 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-510',
    },
    'raw_irradiance_555': {
        'long_name': 'Raw Downwelling Irradiance 555 nm',
        'comment': ('Raw downwelling spectral irradiance at 555 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-555',
    },
    'raw_irradiance_620': {
        'long_name': 'Raw Downwelling Irradiance 620 nm',
        'comment': ('Raw downwelling spectral irradiance at 620 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-620',
    },
    'raw_irradiance_683': {
        'long_name': 'Raw Downwelling Irradiance 683 nm',
        'comment': ('Raw downwelling spectral irradiance at 683 nm measured by the Sea-Bird Electronics (formerly '
                    'Satlantic) OCR-507 sensor.'),
        'units': 'count',
        'data_product_identifier': 'SPECTIR_L0-683',
    },
    # dataset attributes --> derived values
    'downwelling_irradiance_412': {
        'long_name': 'Downwelling Spectral Irradiance at 412 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 412,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 412 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-412',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_412'
    },
    'downwelling_irradiance_444': {
        'long_name': 'Downwelling Spectral Irradiance at 444 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 444,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 444 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-444',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_444'
    },
    'downwelling_irradiance_490': {
        'long_name': 'Downwelling Spectral Irradiance at 490 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 490,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 490 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-490',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_490'
    },
    'downwelling_irradiance_510': {
        'long_name': 'Downwelling Spectral Irradiance at 510 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 510,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 510 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-510',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_510'
    },
    'downwelling_irradiance_555': {
        'long_name': 'Downwelling Spectral Irradiance at 555 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 555,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 555 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-555',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_555'
    },
    'downwelling_irradiance_620': {
        'long_name': 'Downwelling Spectral Irradiance at 620 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 620,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 620 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-620',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_620'
    },
    'downwelling_irradiance_683': {
        'long_name': 'Downwelling Spectral Irradiance at 683 nm',
        'standard_name': 'downwelling_photon_spherical_irradiance_per_unit_wavelength_in_sea_water',
        'radiation_wavelength': 683,
        'units': 'uW cm-2 nm-1',
        'comment': ('Downwelling spectral irradiance measured at 683 nm by the Sea-Bird Electronics OCR-507 '
                    'Multispectral Radiometer. Spectral irradiance is a critical measurement for defining important '
                    'ocean processes, such as the radiant heating rate, and sets the energy available to drive a '
                    'range of biological and chemical processes in the ocean.'),
        'data_product_identifier': 'SPECTIR_L1-683',
        '_FillValue': np.nan,
        'ancillary_variables': 'raw_irradiance_683'
    }
})


def spkir_datalogger(ds, burst=False):
    """
    Takes spkir data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly. Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial spkir data set downloaded from OOI via the M2M system
    :param burst: resample the data to the defined time interval
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   instrument_id == captured in global metadata (always the OCR-507)
    #   serial_number == captured in global metadata
    #   dcl_controller_timestamp == empty string
    #   passed_checksum == not ingested if it doesn't pass, don't need to restate that here
    #   internal_timestamp == superseded by time, redundant so can remove
    #   frame_counter = not used
    #   sample_delay == superseded by time, redundant so can remove
    #   timer == superseded by time, redundant so can remove
    drop_vars = ['instrument_id', 'serial_number', 'dcl_controller_timestamp', 'passed_checksum',
                 'internal_timestamp', 'frame_counter', 'sample_delay', 'timer']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # rename some parameters here to get a better defined data set with cleaner attributes
    rename = {
        'vin_sense': 'input_voltage',
        'va_sense': 'analog_rail_voltage'
    }
    ds = ds.rename(rename)

    # convert raw voltages and the instrument temperature to engineering units
    ds['input_voltage'] = ds['input_voltage'] * 0.03
    ds['analog_rail_voltage'] = ds['analog_rail_voltage'] * 0.03
    ds['internal_temperature'] = -50 + ds['internal_temperature'] * 0.5

    # add the raw downwelling spectral irradiance measurements to the data set split out by wavelength
    raw = ds['channel_array']
    ds['raw_irradiance_412'] = raw[:, 0]
    ds['raw_irradiance_444'] = raw[:, 1]
    ds['raw_irradiance_490'] = raw[:, 2]
    ds['raw_irradiance_510'] = raw[:, 3]
    ds['raw_irradiance_555'] = raw[:, 4]
    ds['raw_irradiance_620'] = raw[:, 5]
    ds['raw_irradiance_683'] = raw[:, 6]

    # add the converted downwelling spectral irradiance measurements to the data set split out by wavelength
    ed = ds['spkir_abj_cspp_downwelling_vector']
    ds['downwelling_irradiance_412'] = ed[:, 0]
    ds['downwelling_irradiance_444'] = ed[:, 1]
    ds['downwelling_irradiance_490'] = ed[:, 2]
    ds['downwelling_irradiance_510'] = ed[:, 3]
    ds['downwelling_irradiance_555'] = ed[:, 4]
    ds['downwelling_irradiance_620'] = ed[:, 5]
    ds['downwelling_irradiance_683'] = ed[:, 6]

    # remove the original 2D variables
    ds = ds.drop_vars(['spectra', 'channel_array', 'spkir_abj_cspp_downwelling_vector'])

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # add the original variable name as an attribute, if renamed
    for key, value in rename.items():
        ds[value].attrs['ooinet_variable_name'] = key

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

    if burst:   # re-sample the data to a defined time interval using a median average
        # create the burst averaging
        ds['time'] = ds['time'] + np.timedelta64(450, 's')
        burst = ds.resample(time='900s', skipna=True).median(dim='time', keep_attrs=True)
        burst = burst.where(~np.isnan(burst.deployment), drop=True)

        # reset the attributes...which keep_attrs should do...
        burst.attrs = ds.attrs
        for v in burst.variables:
            burst[v].attrs = ds[v].attrs

        # save the newly average data
        ds = burst

    return ds


def spkir_cspp(ds):
    """
    Takes SPKIR data recorded by the CSPP loggers used by the Endurance Array
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial SPKIR data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   instrument_id == captured in global metadata (always the OCR-507)
    #   internal_timestamp == time, redundant so can remove
    #   profiler_timestamp == time, redundant so can remove
    #   suspect_timestamp = not used
    #   frame_counter = not used
    #   sample_delay == superseded by time, redundant so can remove
    #   timer == superseded by time, redundant so can remove
    drop_vars = ['instrument_id', 'internal_timestamp', 'profiler_timestamp',
                 'suspect_timestamp', 'frame_counter', 'sample_delay', 'timer']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # rename here for consistency across other data sets
    rename = {
        'vin_sense': 'input_voltage',
        'va_sense': 'analog_rail_voltage',
        'pressure': 'seawater_pressure',
        'pressure_qc_executed': 'seawater_pressure_qc_executed',
        'pressure_qc_results': 'seawater_pressure_qc_results',
    }
    ds = ds.rename(rename)

    # convert raw voltages and the instrument temperature to engineering units
    ds['input_voltage'] = ds['input_voltage'] * 0.03
    ds['analog_rail_voltage'] = ds['analog_rail_voltage'] * 0.03
    ds['internal_temperature'] = -50 + ds['internal_temperature'] * 0.5

    # add the raw downwelling spectral irradiance measurements to the data set split out by wavelength
    raw = ds['channel_array']
    ds['raw_irradiance_412'] = raw[:, 0]
    ds['raw_irradiance_444'] = raw[:, 1]
    ds['raw_irradiance_490'] = raw[:, 2]
    ds['raw_irradiance_510'] = raw[:, 3]
    ds['raw_irradiance_555'] = raw[:, 4]
    ds['raw_irradiance_620'] = raw[:, 5]
    ds['raw_irradiance_683'] = raw[:, 6]

    # add the converted downwelling spectral irradiance measurements to the data set split out by wavelength
    ed = ds['spkir_abj_cspp_downwelling_vector']
    ds['downwelling_irradiance_412'] = ed[:, 0]
    ds['downwelling_irradiance_444'] = ed[:, 1]
    ds['downwelling_irradiance_490'] = ed[:, 2]
    ds['downwelling_irradiance_510'] = ed[:, 3]
    ds['downwelling_irradiance_555'] = ed[:, 4]
    ds['downwelling_irradiance_620'] = ed[:, 5]
    ds['downwelling_irradiance_683'] = ed[:, 6]

    # remove the original 2D variables
    ds = ds.drop_vars(['spectra', 'channel_array', 'spkir_abj_cspp_downwelling_vector'])

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # add the original variable name as an attribute, if renamed
    for key, value in rename.items():
        ds[value].attrs['ooinet_variable_name'] = key

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)

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
    burst = args.burst

    # check if we are specifying a deployment or a specific date and time range
    if not deploy or (start and stop):
        return SyntaxError('You must specify either a deployment number or beginning and end dates of interest.')

    # if we are specifying a deployment number, then get the data from the Gold Copy THREDDS server
    if deploy:
        # download the data for the deployment
        spkir = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*SPKIR.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not spkir:
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
        spkir = m2m_collect(r, '.*SPKIR.*\\.nc$')

        # check to see if we downloaded any data
        if not spkir:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize the data
    if node == 'SP001':
        # this SPKIR is part of a CSPP
        spkir = spkir_cspp(spkir)
    else:
        # this SPKIR is standalone on one of the Surface Moorings
        spkir = spkir_datalogger(spkir, burst)

    vocab = get_vocabulary(site, node, sensor)[0]
    spkir = update_dataset(spkir, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    spkir.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
