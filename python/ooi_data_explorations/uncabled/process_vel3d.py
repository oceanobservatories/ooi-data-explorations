#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os

from ooi_data_explorations.common import inputs, load_gc_thredds, m2m_collect, m2m_request, \
    update_dataset, ENCODINGS
from ooi_data_explorations.qartod.qc_processing import parse_qc

# load configuration settings
FILL_INT = -9999999
ATTRS = dict({
    'noise_amplitude_beam1': {
        'long_name': 'Noise Amplitude Beam 1',
        'comment': ('Ambient noise amplitudes measured by beam 1 prior to a the collection of a 3 minute burst of '
                    '8 Hz velocity data and recorded in the header data packet. For the vector, these values '
                    'should be around 50.'),
        'units': 'counts'
    },
    'noise_amplitude_beam2': {
        'long_name': 'Noise Amplitude Beam 2',
        'comment': ('Ambient noise amplitudes measured by beam 2 prior to a the collection of a 3 minute burst of '
                    '8 Hz velocity data and recorded in the header data packet. For the vector, these values '
                    'should be around 50.'),
        'units': 'counts'
    },
    'noise_amplitude_beam3': {
        'long_name': 'Noise Amplitude Beam 3',
        'comment': ('Ambient noise amplitudes measured by beam 3 prior to a the collection of a 3 minute burst of '
                    '8 Hz velocity data and recorded in the header data packet. For the vector, these values '
                    'should be around 50.'),
        'units': 'counts'
    },
    'noise_correlation_beam1': {
        'long_name': 'Noise Correlation Beam 1',
        'comment': ('Ambient noise correlations measured by beam 1 prior to the collection of a 3 minute burst of '
                    '8 Hz velocity data and recorded in the header data packet. The correlation is a measure of the '
                    'similarity of consecutive measurements. Zero correlation means nothing at all is similar between '
                    'the measurements (e.g. noise is all over the place), whereas a correlation of 100 means the '
                    'measurements are identical. We want high correlations because it gives us confidence the system '
                    'measured the ambient noise field correctly.'),
        'units': 'percent'
    },
    'noise_correlation_beam2': {
        'long_name': 'Noise Correlation Beam 2',
        'comment': ('Ambient noise correlations measured by beam 2 prior to the collection of a 3 minute burst of '
                    '8 Hz velocity data and recorded in the header data packet. The correlation is a measure of the '
                    'similarity of consecutive measurements. Zero correlation means nothing at all is similar between '
                    'the measurements (e.g. noise is all over the place), whereas a correlation of 100 means the '
                    'measurements are identical. We want high correlations because it gives us confidence the system '
                    'measured the ambient noise field correctly.'),
        'units': 'percent'
    },
    'noise_correlation_beam3': {
        'long_name': 'Noise Correlation Beam 3',
        'comment': ('Ambient noise correlations measured by beam 3 prior to the collection of a 3 minute burst of '
                    '8 Hz velocity data and recorded in the header data packet. The correlation is a measure of the '
                    'similarity of consecutive measurements. Zero correlation means nothing at all is similar between '
                    'the measurements (e.g. noise is all over the place), whereas a correlation of 100 means the '
                    'measurements are identical. We want high correlations because it gives us confidence the system '
                    'measured the ambient noise field correctly.'),
        'units': 'percent'
    },
    'battery_voltage': {
        'long_name': 'Battery Voltage',
        'comment': 'Voltage of either the internal battery pack or externally supplied power, whichever is greater.',
        'units': 'V'
    },
    'speed_of_sound': {
        'long_name': 'Speed of Sound',
        'comment': ('Estimated speed of sound derived internally by the VEL3D from the temperature sensor '
                    'measurements and an assumed constant salinity of 33 psu.'),
        'units': 'm s-1'
    },
    'heading': {
        'long_name': 'Heading',
        'comment': 'Measured heading of the VEL3D, uncorrected for magnetic declination.',
        'units': 'degrees'
    },
    'pitch': {
        'long_name': 'Pitch',
        'comment': 'Measured pitch of the VEL3D.',
        'units': 'degrees'
    },
    'roll': {
        'long_name': 'Roll',
        'comment': 'Measured roll of the VEL3D.',
        'units': 'degrees'
    },
    'temperature': {
        'long_name': 'Sea Water Temperature',
        'standard_name': 'sea_water_temperature',
        'comment': 'In-situ sea water temperature measured at the base of the transducer stalk.',
        'units': 'degrees_Celsius'
    },
    'error_code': {
        'long_name': 'Instrument Error Codes',
        'flag_masks': np.array([1, 2, 4, 8, 16, 32, 64], dtype=np.uint8),
        'flag_meanings': ('compass_error measurement_error sensor_data_error tag_bit_error '
                          'flash_error undefined ct_sensor_read_error'),
        'comment': 'Integer representation of the instrument error codes.'
        # 'units': '',    deliberately left blank, no units for this value
    },
    'status_code': {
        'long_name': 'Instrument Status Codes',
        'flag_masks': np.array([1, 2, 4, 8, 48, 48, 48, 48, 192, 192, 192, 192], dtype=np.uint8),
        'flag_values': np.array([1, 2, 4, 8, 0, 16, 32, 48, 0, 64, 128, 192], dtype=np.uint8),
        'flag_meanings': ('orientation_down scaling_factor_0.1 pitch_out_of_range roll_out_of_range '
                          'wake_bad_power wake_break_received wake_power_applied wake_rtc_alarm '
                          'power_level_high power_level_1 power_level_2 power_level_low'),
        'comment': 'Integer representation of the instrument status codes.'
        # 'units': '',    deliberately left blank, no units for this value
    },
    'pressure': {
        'long_name': 'Pressure',
        'standard_name': 'sea_water_pressure_due_to_sea_water',
        'comment': 'Sea water pressure measured at the base of the transducer stalk.',
        'units': 'dbar'
    },
    'velocity_east': {
        'long_name': 'Estimated Eastward Sea Water Velocity',
        'comment': 'Estimated eastward sea water velocity uncorrected for magnetic declination.',
        'data_product_identifier': 'VELPTTU-VLE_L0',
        'units': 'mm s-1',
    },
    'velocity_east_corrected': {
        'long_name': 'Eastward Sea Water Velocity',
        'standard_name': 'eastward_sea_water_velocity',
        'comment': 'Eastward sea water velocity corrected for magnetic declination and scaled to m/s.',
        'data_product_identifier': 'VELPTTU-VLE_L1',
        'units': 'm s-1',
    },
    'velocity_north': {
        'long_name': 'Estimated Northward Sea Water Velocity',
        'comment': 'Estimated northward sea water velocity uncorrected for magnetic declination.',
        'data_product_identifier': 'VELPTTU-VLN_L0',
        'units': 'mm s-1',
    },
    'velocity_north_corrected': {
        'long_name': 'Northward Sea Water Velocity',
        'standard_name': 'northward_sea_water_velocity',
        'comment': 'Northward sea water velocity corrected for magnetic declination and scaled to m/s.',
        'data_product_identifier': 'VELPTTU-VLN_L1',
        'units': 'm s-1',
    },
    'velocity_vertical': {
        'long_name': 'Upward Sea Water Velocity',
        'standard_name': 'upward_sea_water_velocity',
        'comment': 'Vertical sea water velocity component.',
        'data_product_identifier': 'VELPTTU-VLU_L0',
        'units': 'mm s-1',
    },
    'amplitude_beam1': {
        'long_name': 'Velocity Amplitude Beam 1',
        'comment': ('Raw measurement, for beam 1, of the difference in frequency between the transmitted '
                    'and the received pulse, which is proportional to the velocity of the water.'),
        'units': 'count'
    },
    'amplitude_beam2': {
        'long_name': 'Velocity Amplitude Beam 2',
        'comment': ('Raw measurement, for beam 2, of the difference in frequency between the transmitted '
                    'and the received pulse, which is proportional to the velocity of the water.'),
        'units': 'count'
    },
    'amplitude_beam3': {
        'long_name': 'Velocity Amplitude Beam 3',
        'comment': ('Raw measurement, for beam 3, of the difference in frequency between the transmitted '
                    'and the received pulse, which is proportional to the velocity of the water.'),
        'units': 'count'
    },
    'correlation_beam1': {
        'long_name': 'Percent Correlation Beam 1',
        'comment': ('Percent correlation, for beam 1, is a measure of the similarity of two pulse echoes being '
                    'measured by the Doppler instrument. Zero correlation means nothing at all is similar between '
                    'the two echoes, whereas a correlation of 100 means the two echoes are identical. We want high '
                    'correlation because it gives us confidence the system measured the two pulses it originally '
                    'sent out and is determining a valid phase shift.'),
        'units': 'percent'
    },
    'correlation_beam2': {
        'long_name': 'Percent Correlation Beam 2',
        'comment': ('Percent correlation, for beam 2, is a measure of the similarity of two pulse echoes being '
                    'measured by the Doppler instrument. Zero correlation means nothing at all is similar between '
                    'the two echoes, whereas a correlation of 100 means the two echoes are identical. We want high '
                    'correlation because it gives us confidence the system measured the two pulses it originally '
                    'sent out and is determining a valid phase shift.'),
        'units': 'percent'
    },
    'correlation_beam3': {
        'long_name': 'Percent Correlation Beam 3',
        'comment': ('Percent correlation, for beam 3, is a measure of the similarity of two pulse echoes being '
                    'measured by the Doppler instrument. Zero correlation means nothing at all is similar between '
                    'the two echoes, whereas a correlation of 100 means the two echoes are identical. We want high '
                    'correlation because it gives us confidence the system measured the two pulses it originally '
                    'sent out and is determining a valid phase shift.'),
        'units': 'percent'
    }
})


def quality_checks(ds):
    """
    Quality assessment of the pitch, roll, and pressure values for the VEL3D
    using a susbset of the QARTOD flags to indicate the quality. QARTOD
    flags used are:

        1 = Pass
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

    # test for valid speed of sound values (between 1400 and 1700 m/s).  This is a very rough test, but should
    # catch the most egregious errors. The speed of sound is calculated from the temperature and nominal salinity
    # values, so this test is really just a sanity check.
    m = ds['speed_of_sound'] <= 1400 or ds['speed_of_sound'] >= 1700
    qc_flag[m] = 4

    # test for pressure out of range (catch those periods when the instrument is out of the water)
    m = ds['seawater_pressure'] <= 20
    qc_flag[m] = 4

    # test for periods when the measured amplitudes are too low compared to the noise amplitudes
    m = ds['amplitude_beam1'] < ds['noise_amplitude_beam1'] * 0.5
    qc_flag[m] = 4

    return qc_flag


def vel3d_datalogger(header, system, velocity, burst=False):
    """
    Takes VEL3D (Nortek Vector) data recorded by the data loggers used in the
    EA moorings and cleans up the data set to make it more user-friendly.
    Primary task is adding the header and system packet data to the velocity
    data set in order to make a single data set with all the information
    needed to use and assess the data. Additionally, this function is used to
    rename parameters and drop some that are of limited use.

    :param header: xarray dataset with the header data. This data set contains
        the ambient noise measurements made prior to the instrument collecting
        a burst of data.
    :param system: xarray dataset with the system data. This data set contains
        the 1 Hz instrument status data (e.g. battery voltage, temperature,
        pressure, tilt, roll, etc.).
    :param velocity: xarray dataset with the velocity data. This data set
        contains the 8 Hz velocity data collected by the instrument.
    :param burst: boolean, used to indicate if the data should be burst
        averaged or not. Default is False.
    :return ds: cleaned up data set
    """
    # drop some of the variables: velocity
    #   analog_input_1 == not used, no data
    #   analog_input_2 == not used, no data
    #   sea_water_pressure_mbar_qc_executed == drop in favor of instrument specific tests
    #   sea_water_pressure_mbar_qc_results == drop in favor of instrument specific tests
    #   turbulent_velocity_east_qc_executed == drop in favor of instrument specific tests
    #   turbulent_velocity_east_qc_results == drop in favor of instrument specific tests
    #   turbulent_velocity_north_qc_executed == drop in favor of instrument specific tests
    #   turbulent_velocity_north_qc_results == drop in favor of instrument specific tests
    #   turbulent_velocity_vertical_qc_executed == drop in favor of instrument specific tests
    #   turbulent_velocity_vertical_qc_results == drop in favor of instrument specific tests
    #   vel3d_c_eastward_turbulent_velocity_qc_executed == drop in favor of instrument specific tests
    #   vel3d_c_eastward_turbulent_velocity_qc_results == drop in favor of instrument specific tests
    #   vel3d_c_northward_turbulent_velocity_qc_executed == drop in favor of instrument specific tests
    #   vel3d_c_northward_turbulent_velocity_qc_results == drop in favor of instrument specific tests
    #   vel3d_c_turbulent_eastward_velocity == velocity_vertical == redundant, so can remove
    #   vel3d_c_upward_turbulent_velocity_qc_executed == drop in favor of instrument specific tests
    #   vel3d_c_upward_turbulent_velocity_qc_results == drop in favor of instrument specific tests
    drop_vars = ['analog_input_1', 'analog_input_2', 'sea_water_pressure_mbar_qc_executed',
                 'sea_water_pressure_mbar_qc_results', 'turbulent_velocity_east_qc_executed',
                 'turbulent_velocity_east_qc_results', 'turbulent_velocity_north_qc_executed',
                 'turbulent_velocity_north_qc_results', 'turbulent_velocity_vertical_qc_executed',
                 'turbulent_velocity_vertical_qc_results', 'vel3d_c_eastward_turbulent_velocity_qc_executed',
                 'vel3d_c_eastward_turbulent_velocity_qc_results', 'vel3d_c_northward_turbulent_velocity_qc_executed',
                 'vel3d_c_northward_turbulent_velocity_qc_results', 'vel3d_c_upward_turbulent_velocity',
                 'vel3d_c_upward_turbulent_velocity_qc_executed', 'vel3d_c_upward_turbulent_velocity_qc_results'
                 ]
    for var in velocity.variables:
        if var in drop_vars:
            velocity = velocity.drop_vars(var)

    # drop some of the variables: system
    #   analog_input == not used, no data
    #   date_time_string == redundant, so can remove
    #   speed_of_sound_qc_executed == drop in favor of instrument specific tests
    #   speed_of_sound_qc_results == drop in favor of instrument specific tests
    drop_vars = ['analog_input', 'date_time_string', 'speed_of_sound_qc_executed', 'speed_of_sound_qc_results']
    for var in system.variables:
        if var in drop_vars:
            system = system.drop_vars(var)

    # drop some of the variables: header
    #   date_time_string == redundant, so can remove
    #   number_velocity_records == this is always 1440, so not useful
    drop_vars = ['date_time_string', 'number_velocity_records']
    for var in header.variables:
        if var in drop_vars:
            header = header.drop_vars(var)

    # rename some parameters here to get a better defined data set with cleaner attributes
    rename = {
        # header packets
        'noise_amp_beam1': 'noise_amplitude_beam1',
        'noise_amp_beam2': 'noise_amplitude_beam2',
        'noise_amp_beam3': 'noise_amplitude_beam3',
        # system packets
        'battery_voltage_dv': 'battery_voltage',
        'heading_decidegree': 'heading',
        'pitch_decidegree': 'pitch',
        'roll_decidegree': 'roll',
        'temperature_centidegree': 'seawater_temperature',
        # velocity packets
        'sea_water_pressure_mbar': 'seawater_pressure',
        'amplitude_beam_1': 'amplitude_beam1',
        'amplitude_beam_2': 'amplitude_beam2',
        'amplitude_beam_3': 'amplitude_beam3',
        'correlation_beam_1': 'correlation_beam1',
        'correlation_beam_2': 'correlation_beam2',
        'correlation_beam_3': 'correlation_beam3',
        'turbulent_velocity_east': 'velocity_east',
        'turbulent_velocity_north': 'velocity_north',
        'turbulent_velocity_vertical': 'velocity_vertical',
        'vel3d_c_turbulent_eastward_velocity': 'corrected_velocity_east',
        'vel3d_c_turbulent_eastward_north': 'corrected_velocity_north',
    }
    for key in rename.keys():
        if key in velocity.variables:
            velocity = velocity.rename({key: rename.get(key)})
        if key in system.variables:
            system = system.rename({key: rename.get(key)})
        if key in header.variables:
            header = header.rename({key: rename.get(key)})

    # convert some variables to more standard units
    system['heading'] = system['heading'] / 10.0  # convert from ddeg to deg
    system['pitch'] = system['pitch'] / 10.0  # convert from ddeg to deg
    system['roll'] = system['roll'] / 10.0  # convert from ddeg to deg
    system['battery_voltage'] = system['battery_voltage'] / 10.0  # convert from dV to V
    system['seawater_temperature'] = system['seawater_temperature'] / 100.0  # convert from cdeg C to deg C

    velocity['seawater_pressure'] = velocity['seawater_pressure'] / 1000.0  # parser doesn't include this needed scaling term

    system['speed_of_sound'] = system['speed_of_sound'] / 10.0  # convert from dm/s to m/s

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
    ds['vector_sensor_quality_flag'] = quality_checks(ds)
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
        vel3d = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*VEL3D.*\\.nc$' % deploy))

        # check to see if we downloaded any data
        if not vel3d:
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
        vel3d = m2m_collect(r, '.*VEL3D.*\\.nc$')

        # check to see if we downloaded any data
        if not vel3d:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # clean-up and reorganize the data
    if node == 'SP001':
        # this VEL3D is part of a CSPP
        vel3d = vel3d_cspp(vel3d)
    else:
        # this VEL3D is standalone on one of the Surface Moorings (one of those rare cases where all forms of the data
        # are in the same format)
        vel3d = vel3d_datalogger(vel3d)

    depth = vel3d.depth.mean(dim='time').values
    vel3d = update_dataset(vel3d, depth)

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    vel3d.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
