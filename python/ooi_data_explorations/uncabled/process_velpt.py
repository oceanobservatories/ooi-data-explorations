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
    'error_code': {
        'long_name': 'Error Code',
        # 'units': '',    deliberately left blank, no units for this value
        'comment': 'Instrument error codes.',
        'flag_mask': np.array(2 ** np.array(range(0, 8)), dtype=object).astype(np.intc),
        'flag_meanings': ('compass_error measurement_data_error sensor_data_error tag_bit_error flash_error '
                          'undefined serial_ct_sensor_read_error undefined')
    },
    'battery_voltage': {
        'long_name': 'Battery Voltage',
        'units': 'V',
        'comment': ('Reports either the internal battery voltage, or the external power applied, whichever is '
                    'greater. Recorded in dV and converted to V.')
    },
    'speed_of_sound': {
        'long_name': 'Speed of Sound',
        'units': 'm s-1',
        'comment': 'Contains either manual or calculated speed of sound recorded in dm/s and converted to m/s.'
    },
    'heading': {
        'long_name': 'Heading',
        'units': 'degrees',
        'comment': ('Measured heading of the instrument, uncorrected for magnetic declination, recorded in '
                    'ddegrees and converted to degrees.')
    },
    'pitch': {
        'long_name': 'Pitch',
        'units': 'degrees',
        'comment': 'Measured pitch of the instrument recorded in ddegrees and converted to degrees.'
    },
    'roll': {
        'long_name': 'Roll',
        'units': 'degrees',
        'comment': 'Measured roll of the instrument recorded in ddegrees and converted to degrees.'
    },
    'seawater_pressure': {
        'long_name': 'Seawater Pressure',
        'standard_name': 'sea_water_pressure_due_to_sea_water',
        'units': 'dbar',
        'comment': 'Instrument pressure sensor value recorded in 0.1 mbar and converted to dbar.',
    },
    'status_code': {
        'long_name': 'Status Code',
        # 'units': '',    deliberately left blank, no units for this value
        'comment': 'Instrument status codes.',
        'flag_mask': np.array([1, 1, 2, 2, 4, 8, 48, 48, 48, 48, 192, 192, 192, 192], dtype=object).astype(np.intc),
        'flag_values': np.array([0, 1, 0, 2, 4, 8, 0, 16, 32, 48, 0, 64, 128, 192], dtype=object).astype(np.intc),
        'flag_meanings': ('orientation_up orientation_down not_used not_used pitch_out_of_range roll_out_of_range '
                          'wakeup_state_bad_power wakeup_state_break wakeup_state_power_applied wakeup_state_rtc_alarm '
                          'power_level_high power_level_mid_high power_level_mid_low power_level_low')
        # Per https://cfconventions.org/cf-conventions/cf-conventions.html#flags, a blend of flag masks, flag values,
        # and flag meanings are used "to describe a blend of independent Boolean conditions and enumerated status
        # codes."
    },
    'seawater_temperature': {
        'long_name': 'Sea Water Temperature',
        'standard_name': 'sea_water_temperature',
        'units': 'degrees_Celsius',
        'comment': ('In-situ sea water temperature measured at the transducer face. Recorded in cdegrees Celsius and '
                    'converted to degrees Celsius.')
    },
    'velocity_east': {
        'long_name': 'Estimated Eastward Seawater Velocity',
        'units': 'mm s-1',
        'comment': ('This is the eastward seawater velocity component uncorrected for magnetic declination as '
                    'reported by the instrument in mm/s.'),
        'data_product_identifier': 'VELPTMN-VLE_L0'
    },
    'velocity_north': {
        'long_name': 'Estimated Northward Seawater Velocity',
        'units': 'mm s-1',
        'comment': ('This is the northward seawater velocity component uncorrected for magnetic declination as '
                    'reported by the instrument in mm/s.'),
        'data_product_identifier': 'VELPTMN-VLN_L0'
    },
    'velocity_vertical': {
        'long_name': 'Upward Seawater Velocity',
        'standard_name': 'upward_sea_water_velocity',
        'units': 'mm s-1',
        'comment': 'The vertical seawater velocity component as reported by the instrument in mm/s.',
        'data_product_identifier': 'VELPTMN-VLU_L0'
    },
    'amplitude_beam1': {
        'long_name': 'Amplitude Beam 1',
        'units': 'count',
        'comment': ('This is the raw measurement, the acoustic return signal for the beam, used to calculate the '
                    'seawater velocity. This value should be roughly equivalent to the other two beams. Significant '
                    'differences would suggest one or more of the beams is blocked or otherwise impaired.')
    },
    'amplitude_beam2': {
        'long_name': 'Amplitude Beam 2',
        'units': 'count',
        'comment': ('This is the raw measurement, the acoustic return signal for the beam, used to calculate the '
                    'seawater velocity. This value should be roughly equivalent to the other two beams. Significant '
                    'differences would suggest one or more of the beams is blocked or otherwise impaired.')
    },
    'amplitude_beam3': {
        'long_name': 'Amplitude Beam 3',
        'units': 'count',
        'comment': ('This is the raw measurement, the acoustic return signal for the beam, used to calculate the '
                    'seawater velocity. This value should be roughly equivalent to the other two beams. Significant '
                    'differences would suggest one or more of the beams is blocked or otherwise impaired.')
    },
    # ---- derived values ----
    'depth': {
        'long_name': 'Depth',
        'standard_name': 'depth',
        'units': 'm',
        'comment': ('Depth of the instrument below the surface of the water. Calculated from the instrument pressure'
                    'sensor.'),
        'positive': 'down'
    },
    'eastward_seawater_velocity': {
        'long_name': 'Eastward Seawater Velocity',
        'standard_name': 'eastward_sea_water_velocity',
        'units': 'm s-1',
        'comment': ('Eastward sea water velocity component in Earth coordinates corrected for magnetic declination '
                    'and scaled to standard units of m s-1.'),
        'data_product_identifier': 'VELPTMN-VLE_L1',
        'ancillary_variables': 'velocity_east, velocity_north, time, lat, lon, z'
    },
    'northward_seawater_velocity': {
        'long_name': 'Northward Seawater Velocity',
        'standard_name': 'northward_sea_water_velocity',
        'units': 'm s-1',
        'comment': ('Northward sea water velocity component in Earth coordinates corrected for magnetic declination '
                    'and scaled to standard units of m s-1.'),
        'data_product_identifier': 'VELPTMN-VLN_L1',
        'ancillary_variables': 'velocity_east, velocity_north, time, lat, lon, z'
    },
    'upward_seawater_velocity': {
        'long_name': 'Upward Seawater Velocity',
        'standard_name': 'upward_sea_water_velocity',
        'units': 'm s-1',
        'comment': ('The vertical seawater velocity component as reported by the instrument in mm/s and scaled to '
                    'm s-1.'),
        'data_product_identifier': 'VELPTMN-VLU_L1',
        'ancillary_variables': 'velocity_vertical'
    }
})


def quality_checks(ds):
    """
    Quality assessment of the pitch, roll, and pressure values for the VELPT
    using a susbset of the QARTOD flags to indicate the quality. QARTOD
    flags used are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail

    The final flag value represents the worst case assessment of the data quality.

    :param ds: xarray dataset with the pitch, roll, and seawater pressure
    :return qc_flag: array of flag values indicating data quality
    """
    qc_flag = ds['time'].astype('int32') * 0 + 1   # default flag values, no errors

    # test for pitch and roll out of range (greater than 30 degrees)
    m = np.abs(ds['pitch']) > 20
    qc_flag[m] = 3
    m = np.abs(ds['pitch']) >= 30
    qc_flag[m] = 4

    m = np.abs(ds['roll']) > 20
    qc_flag[m] = 3
    m = np.abs(ds['roll']) >= 30
    qc_flag[m] = 4

    # test for pressure out of range (catch those periods when the instrument is out of the water)
    if ds.attrs['node'] in ['SBD11', 'SBD12', 'SBD17']:
        # surface buoy, pressure is nominally 1.25 dbar
        m = ds['seawater_pressure'] <= 0
    else:
        # subsurface platform (NSIF or MFN), pressure is always greater than 3 dbar
        m = ds['seawater_pressure'] <= 3

    qc_flag[m] = 4

    return qc_flag


def velpt_datalogger(ds):
    """
    Takes VELPT (Nortek Aquadopp) data recorded by the data loggers used in the
    CGSN/EA moorings and cleans up the data set to make it more user-friendly.
    Primary task is renaming parameters and dropping some that are of limited
    use. Additionally, re-organize some of the variables to permit better
    assessments of the data.

    :param ds: initial velpt data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   analog1 == not used, not data
    #   internal_timestamp == time, redundant so can remove
    #   date_time_string == time, redundant so can remove
    #   velocity_beam1_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam1_qc_results == QC tests are not applied to L0 data
    #   velocity_beam2_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam2_qc_results == QC tests are not applied to L0 data
    #   velocity_beam3_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam3_qc_results == QC tests are not applied to L0 data
    drop_vars = ['analog1', 'internal_timestamp', 'date_time_string',
                 'velocity_beam1_qc_executed', 'velocity_beam1_qc_results',
                 'velocity_beam2_qc_executed', 'velocity_beam2_qc_results',
                 'velocity_beam3_qc_executed', 'velocity_beam3_qc_results'
                 ]
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # rename some parameters here to get a better defined data set with cleaner attributes
    rename = {
        'battery_voltage_dv': 'battery_voltage',
        'heading_decidegree': 'heading',
        'pitch_decidegree': 'pitch',
        'roll_decidegree': 'roll',
        'status': 'status_code',
        'sea_water_pressure_mbar': 'seawater_pressure',
        'sea_water_pressure_mbar_qc_executed': 'seawater_pressure_qc_executed',
        'sea_water_pressure_mbar_qc_results': 'seawater_pressure_qc_results',
        'sound_speed_dms': 'speed_of_sound',
        'temperature_centidegree': 'seawater_temperature',
        'velocity_beam1': 'velocity_east',
        'velocity_beam2': 'velocity_north',
        'velocity_beam3': 'velocity_upward',
        'eastward_velocity': 'eastward_seawater_velocity',
        'eastward_velocity_qc_executed': 'eastward_seawater_velocity_qc_executed',
        'eastward_velocity_qc_results': 'eastward_seawater_velocity_qc_results',
        'northward_velocity': 'northward_seawater_velocity',
        'northward_velocity_qc_executed': 'northward_seawater_velocity_qc_executed',
        'northward_velocity_qc_results': 'northward_seawater_velocity_qc_results',
        'upward_velocity': 'upward_seawater_velocity',
        'upward_velocity_qc_executed': 'upward_seawater_velocity_qc_executed',
        'upward_velocity_qc_results': 'upward_seawater_velocity_qc_results'
    }
    ds = ds.rename(rename)

    # convert some variables to more standard units
    ds['heading'] = ds['heading'] / 10.0  # convert from ddeg to deg
    ds['pitch'] = ds['pitch'] / 10.0  # convert from ddeg to deg
    ds['roll'] = ds['roll'] / 10.0  # convert from ddeg to deg
    ds['battery_voltage'] = ds['battery_voltage'] / 10.0  # convert from dV to V
    ds['seawater_pressure'] = ds['seawater_pressure'] / 1000.0  # parser doesn't include this needed scaling term
    ds['seawater_temperature'] = ds['seawater_temperature'] / 100.0  # convert from cdeg C to deg C
    ds['speed_of_sound'] = ds['speed_of_sound'] / 10.0  # convert from dm/s to m/s

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

    # test the data quality using additional instrument variables
    ds['aquadopp_sensor_quality_flag'] = quality_checks(ds)
    return ds


def velpt_instrument(ds):
    """
    Takes VELPT (Nortek Aquadopp) data recorded internally by the instrument
    used and cleans up the data set to make it more user-friendly.
    Primary task is renaming parameters and dropping some that are of limited
    use. Additionally, re-organize some of the variables to permit better
    assessments of the data.

    :param ds: initial velpt data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   analog1 == not used, not data
    #   internal_timestamp == time, redundant so can remove
    #   date_time_string == time, redundant so can remove
    #   velocity_beam1_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam1_qc_results == QC tests are not applied to L0 data
    #   velocity_beam2_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam2_qc_results == QC tests are not applied to L0 data
    #   velocity_beam3_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam3_qc_results == QC tests are not applied to L0 data
    drop_vars = ['analog1', 'internal_timestamp', 'date_time_string',
                 'velocity_beam1_qc_executed', 'velocity_beam1_qc_results',
                 'velocity_beam2_qc_executed', 'velocity_beam2_qc_results',
                 'velocity_beam3_qc_executed', 'velocity_beam3_qc_results'
                 ]
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # rename some parameters here to get a better defined data set with cleaner attributes
    rename = {
        'battery_voltage_dv': 'battery_voltage',
        'heading_decidegree': 'heading',
        'pitch_decidegree': 'pitch',
        'roll_decidegree': 'roll',
        'status': 'status_code',
        'sea_water_pressure_mbar': 'seawater_pressure',
        'sea_water_pressure_mbar_qc_executed': 'seawater_pressure_qc_executed',
        'sea_water_pressure_mbar_qc_results': 'seawater_pressure_qc_results',
        'sound_speed_dms': 'speed_of_sound',
        'temperature_centidegree': 'seawater_temperature',
        'velocity_beam1': 'velocity_east',
        'velocity_beam2': 'velocity_north',
        'velocity_beam3': 'velocity_upward',
        'eastward_velocity': 'eastward_seawater_velocity',
        'eastward_velocity_qc_executed': 'eastward_seawater_velocity_qc_executed',
        'eastward_velocity_qc_results': 'eastward_seawater_velocity_qc_results',
        'northward_velocity': 'northward_seawater_velocity',
        'northward_velocity_qc_executed': 'northward_seawater_velocity_qc_executed',
        'northward_velocity_qc_results': 'northward_seawater_velocity_qc_results',
        'upward_velocity': 'upward_seawater_velocity',
        'upward_velocity_qc_executed': 'upward_seawater_velocity_qc_executed',
        'upward_velocity_qc_results': 'upward_seawater_velocity_qc_results'
    }
    ds = ds.rename(rename)

    # convert some variables to more standard units
    ds['heading'] = ds['heading'] / 10.0  # convert from ddeg to deg
    ds['pitch'] = ds['pitch'] / 10.0  # convert from ddeg to deg
    ds['roll'] = ds['roll'] / 10.0  # convert from ddeg to deg
    ds['battery_voltage'] = ds['battery_voltage'] / 10.0  # convert from dV to V
    ds['seawater_pressure'] = ds['seawater_pressure'] / 1000.0  # parser doesn't include this needed scaling term
    ds['seawater_temperature'] = ds['seawater_temperature'] / 100.0  # convert from cdeg C to deg C
    ds['speed_of_sound'] = ds['speed_of_sound'] / 10.0  # convert from dm/s to m/s

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

    # test the data quality using additional instrument variables
    ds['aquadopp_sensor_quality_flag'] = quality_checks(ds)
    return ds


def velpt_cspp(ds):
    """
    Takes VELPT data recorded by the CSPP loggers used by the Endurance Array
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial VELPT data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   analog1 == not used, not data
    #   internal_timestamp == time, redundant so can remove
    #   profiler_timestamp == time, redundant so can remove
    #   suspect_timestamp = not used
    #   velocity_beam1_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam1_qc_results == QC tests are not applied to L0 data
    #   velocity_beam2_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam2_qc_results == QC tests are not applied to L0 data
    #   velocity_beam3_qc_executed == QC tests are not applied to L0 data
    #   velocity_beam3_qc_results == QC tests are not applied to L0 data
    drop_vars = ['analog1', 'internal_timestamp', 'profiler_timestamp', 'suspect_timestamp',
                 'speed_of_sound_qc_executed', 'speed_of_sound_qc_results',
                 'velocity_beam1_m_s_qc_executed', 'velocity_beam1_m_s_qc_results',
                 'velocity_beam2_m_s_qc_executed', 'velocity_beam2_m_s_qc_results',
                 'velocity_beam3_m_s_qc_executed', 'velocity_beam3_m_s_qc_results',
                 ]
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # rename here for consistency across other data sets
    rename = {
        'velpt_pressure': 'seawater_pressure',
        'velpt_pressure_qc_executed': 'seawater_pressure_qc_executed',
        'velpt_pressure_qc_results': 'seawater_pressure_qc_results',
        'pressure': 'ctd_pressure',
        'pressure_qc_executed': 'ctd_pressure_qc_executed',
        'pressure_qc_results': 'ctd_pressure_qc_results',
        'temperature': 'seawater_temperature',
        'temperature_qc_executed': 'seawater_temperature_qc_executed',
        'temperature_qc_results': 'seawater_temperature_qc_results',
        'velocity_beam1_m_s': 'velocity_east',
        'velocity_beam2_m_s': 'velocity_north',
        'velocity_beam3_m_s': 'velocity_upward',
        'velpt_j_eastward_velocity': 'eastward_seawater_velocity',
        'velpt_j_eastward_velocity_qc_executed': 'eastward_seawater_velocity_qc_executed',
        'velpt_j_eastward_velocity_qc_results': 'eastward_seawater_velocity_qc_results',
        'velpt_j_northward_velocity': 'northward_seawater_velocity',
        'velpt_j_northward_velocity_qc_executed': 'northward_seawater_velocity_qc_executed',
        'velpt_j_northward_velocity_qc_results': 'northward_seawater_velocity_qc_results',
        'velpt_j_upward_velocity': 'upward_seawater_velocity',
        'velpt_j_upward_velocity_qc_executed': 'upward_seawater_velocity_qc_executed',
        'velpt_j_upward_velocity_qc_results': 'upward_seawater_velocity_qc_results',
    }
    ds = ds.rename(rename)

    # correct the velocity variables, which were double-scaled
    ds['velocity_east'] = ds['velocity_east'] * 1000.0  # convert to mm/s from m/s
    ds['velocity_north'] = ds['velocity_north'] * 1000.0  # convert to mm/s from m/s
    ds['velocity_upward'] = ds['velocity_upward'] * 1000.0  # convert to mm/s from m/s
    ds['eastward_seawater_velocity'] = ds['eastward_seawater_velocity'] * 1000.0  # convert to m/s from um/s
    ds['northward_seawater_velocity'] = ds['northward_seawater_velocity'] * 1000.0  # convert to m/s from um/s
    ds['upward_seawater_velocity'] = ds['upward_seawater_velocity'] * 1000.0  # convert to m/s from um/s

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
        velpt = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*VELPT.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not velpt:
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
        velpt = m2m_collect(r, '.*VELPT.*\\.nc$')

        # check to see if we downloaded any data
        if not velpt:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize the data
    if node == 'SP001':
        # this VELPT is part of a CSPP
        velpt = velpt_cspp(velpt)
    else:
        # this VELPT is standalone on one of the Surface Moorings (one of those rare cases where all forms of the data
        # are in the same format)
        velpt = velpt_datalogger(velpt)

    depth = velpt.depth.mean(dim='time').values
    velpt = update_dataset(velpt, depth)

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    velpt.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
