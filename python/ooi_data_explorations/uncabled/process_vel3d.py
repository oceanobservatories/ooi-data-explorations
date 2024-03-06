#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import sys
import re
import xarray as xr

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from tqdm import tqdm

from ooi_data_explorations.common import inputs, get_vocabulary, load_gc_thredds, m2m_collect, m2m_request, \
    get_deployment_dates, update_dataset, ENCODINGS, N_CORES
from ooi_data_explorations.profilers import create_profile_id, bin_profiles
from ooi_data_explorations.qartod.qc_processing import parse_qc

# load configuration settings
FILL_INT = -9999999
VECTOR = dict({
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
                    '8 Hz velocity data and recorded in the header data packet.'),
        'units': 'percent'
    },
    'noise_correlation_beam2': {
        'long_name': 'Noise Correlation Beam 2',
        'comment': ('Ambient noise correlations measured by beam 2 prior to the collection of a 3 minute burst of '
                    '8 Hz velocity data and recorded in the header data packet.'),
        'units': 'percent'
    },
    'noise_correlation_beam3': {
        'long_name': 'Noise Correlation Beam 3',
        'comment': ('Ambient noise correlations measured by beam 3 prior to the collection of a 3 minute burst of '
                    '8 Hz velocity data and recorded in the header data packet.'),
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
    'sea_water_temperature': {
        'long_name': 'Sea Water Temperature',
        'standard_name': 'sea_water_temperature',
        'comment': 'In-situ sea water temperature measured at the base of the transducer stalk.',
        'units': 'degrees_Celsius'
    },
    'error_code': {
        'long_name': 'Instrument Error Codes',
        'flag_masks': np.array([0, 1, 2, 4, 8, 16, 32, 64, 128], dtype=np.uint8),
        'flag_meanings': ('no_errors compass_error measurement_error sensor_data_error tag_bit_error '
                          'flash_error beam_order_error tilt_sensor_error coordinate_transform_error'),
        'comment': 'Integer representation of the instrument error codes.'
        # 'units': '',    deliberately left blank, no units for this value
    },
    'status_code': {
        'long_name': 'Instrument Status Codes',
        'flag_masks': np.array([0, 1, 2, 4, 8, 48, 48, 48, 48, 192, 192, 192, 192], dtype=np.uint8),
        'flag_values': np.array([0, 1, 2, 4, 8, 0, 16, 32, 48, 0, 64, 128, 192], dtype=np.uint8),
        'flag_meanings': ('orientation_up orientation_down scaling_factor_0.1 pitch_out_of_range roll_out_of_range '
                          'wake_bad_power wake_break_received wake_power_applied wake_rtc_alarm power_level_high '
                          'power_level_1 power_level_2 power_level_low'),
        'comment': 'Integer representation of the instrument status codes.'
        # 'units': '',    deliberately left blank, no units for this value
    },
    'scaling_factor': {
        'long_name': 'Scaling Factor',
        'comment': ('Scaling factor used to adjust the velocity data based on bit 1 in the status code, which '
                    'indicates the scaling of the velocity output and depends on the velocity range setting. If the '
                    'instrument is set to use the highest ranges the least significant bit is 1 mm/s. For the lowest '
                    'range it is 0.1 mm/s. This scaling factor has already been applied to the corrected (magnetic '
                    'declination applied) horizontal velocities, but not original the horizontal and vertical '
                    'velocities. Adding it here, so that all velocity measurements are correctly scaled to m/s.'),
        'units': 'mm s-1'
    },
    'sea_water_pressure': {
        'long_name': 'Sea Water Pressure',
        'standard_name': 'sea_water_pressure_due_to_sea_water',
        'comment': 'Sea water pressure measured at the base of the transducer stalk.',
        'units': 'dbar'
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
    },
    'snr_beam1': {
        'long_name': 'Signal-to-Noise Ratio, Beam 1',
        'comment': ('Signal-to-noise ratio, for beam 1, is a measure of the ratio of the signal power to the noise '
                    'power. The higher the ratio, the better the signal. Vector SNR values should be greater than '
                    '15 dB.'),
        'units': 'dB'
    },
    'snr_beam2': {
        'long_name': 'Signal-to-Noise Ratio, Beam 2',
        'comment': ('Signal-to-noise ratio, for beam 2, is a measure of the ratio of the signal power to the noise '
                    'power. The higher the ratio, the better the signal. Vector SNR values should be greater than '
                    '15 dB.'),
        'units': 'dB'
    },
    'snr_beam3': {
        'long_name': 'Signal-to-Noise Ratio, Beam 3',
        'comment': ('Signal-to-noise ratio, for beam 3, is a measure of the ratio of the signal power to the noise '
                    'power. The higher the ratio, the better the signal. Vector SNR values should be greater than '
                    '15 dB.'),
        'units': 'dB'
    },
    'velocity_east': {
        'long_name': 'Estimated Eastward Sea Water Velocity',
        'comment': 'Estimated eastward sea water velocity uncorrected for magnetic declination.',
        'data_product_identifier': 'VELPTTU-VLE_L0',
        'units': 'm s-1',
    },
    'velocity_east_corrected': {
        'long_name': 'Eastward Sea Water Velocity',
        'standard_name': 'eastward_sea_water_velocity',
        'comment': 'Eastward sea water velocity corrected for magnetic declination.',
        'data_product_identifier': 'VELPTTU-VLE_L1',
        'units': 'm s-1',
    },
    'velocity_north': {
        'long_name': 'Estimated Northward Sea Water Velocity',
        'comment': 'Estimated northward sea water velocity uncorrected for magnetic declination.',
        'data_product_identifier': 'VELPTTU-VLN_L0',
        'units': 'm s-1',
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
        'data_product_identifier': 'VELPTTU-VLU_L1',
        'units': 'm s-1',
    },
    'vector_sensor_quality_flag': {
        'long_name': 'Vector Sensor Quality Flag',
        'standard_name': 'aggregate_quality_flag',
        'comment': ('The Vector sensor quality flag is built from a set of vendor defined quality checks looking '
                    'at a combination of different variables to create a single automated quality assessment of the '
                    'data. The final flag value represents the worst case assessment of the data quality. Note, this '
                    'flag represents only an automated quality assessment of the data and should be used in '
                    'conjunction with other flags and data quality information.'),
        'flag_values': np.array([1, 2, 3, 4, 9]),
        'flag_meanings': 'pass not_evaluated suspect_or_of_high_interest fail missing'
    }
})

AQUADOPP = dict({
    'speed_of_sound': {
        'long_name': 'Speed of Sound',
        'comment': ('Estimated speed of sound derived internally by the VEL3D from the temperature sensor '
                    'measurements and an assumed constant salinity of 33 psu.'),
        'units': 'm s-1'
    },
    'sea_water_temperature': {
        'long_name': 'Sea Water Temperature',
        'standard_name': 'sea_water_temperature',
        'comment': 'In-situ sea water temperature measured at the base of the transducer stalk.',
        'units': 'degrees_Celsius'
    },
    'ctd_pressure': {
        'long_name': 'C0-Located CTD Pressure',
        'standard_name': 'sea_water_pressure_due_to_sea_water',
        'comment': ('Sea water pressure from the co-located CTD, interploated into the data record as a more accurate '
                    'pressure sensor.'),
        'units': 'dbar'
    },
    'sea_water_pressure': {
        'long_name': 'Sea Water Pressure',
        'standard_name': 'sea_water_pressure_due_to_sea_water',
        'comment': ('Sea water pressure measured in the center of the transducer face alongside the temperature '
                    'sensor.'),
        'units': 'dbar'
    },
    'heading': {
        'long_name': 'Heading',
        'comment': 'Measured heading of the Aquadopp II, uncorrected for magnetic declination.',
        'units': 'degrees'
    },
    'pitch': {
        'long_name': 'Pitch',
        'comment': 'Measured pitch of the Aquadopp II.',
        'units': 'degrees'
    },
    'roll': {
        'long_name': 'Roll',
        'comment': 'Measured roll of the Aquadopp II.',
        'units': 'degrees'
    },
    'error_code': {
        'long_name': 'Instrument Error Codes',
        'flag_masks': np.array([0, 2, 4, 16, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768], dtype=np.uint16),
        'flag_meanings': ('no_errors data_retrieval_overflow data_retrieval_underrun measurement_error '
                          'tag_error_beam0_in_phase tag_error_beam0_quadrature_phase '
                          'tag_error_beam1_in_phase tag_error_beam1_quadrature_phase '
                          'tag_error_beam2_in_phase tag_error_beam2_quadrature_phase '
                          'tag_error_beam3_in_phase tag_error_beam3_quadrature_phase'),
        'comment': 'Integer representation of the instrument error codes.'
        # 'units': '',    deliberately left blank, no units for this value
    },
    'status_code': {
        'long_name': 'Instrument Status Codes',
        'flag_masks': np.array([48, 48, 48, 48, 192, 192, 192, 192], dtype=np.uint16),
        'flag_values': np.array([0, 16, 32, 48, 0, 64, 128, 192], dtype=np.uint16),
        'flag_meanings': ('wake_bad_power wake_power_applied wake_break_received wake_rtc_alarm '
                          'power_level_high power_level_1 power_level_2 power_level_low'),
        'comment': 'Integer representation of the instrument status codes.'
        # 'units': '',    deliberately left blank, no units for this value
    },
    'cell_size': {
        'long_name': 'Cell Size',
        'comment': 'Length of the cell over which the velocity is measured.',
        'units': 'm'
    },
    'blanking_distance': {
        'long_name': 'Blanking Distance',
        'comment': 'Distance from the transducer face to the first measurement cell.',
        'units': 'm'
    },
    'velocity_range': {
        'long_name': 'Velocity Range',
        'comment': 'Range setting for the velocity measurements. Sets the absolute maximum velocity that can be '
                   'measured. The beam velocity data, once scaled to m/s, can be compared to this value to determine '
                   'if the instrument was operating within its specified range.',
        'units': 'm s-1'
    },
    'battery_voltage': {
        'long_name': 'Battery Voltage',
        'comment': 'Voltage of either the internal battery pack or externally supplied power, whichever is greater.',
        'units': 'V'
    },
    'magnetometer_x': {
        'long_name': 'Magnetometer X',
        'comment': 'Magnetic field strength (Hx) in the x direction. Units are not defined in the manual, but are'
                   'believed to be uT. Combined with the pitch and roll sensor data, the magnetometer data is used'
                   'to determine the heading of the instrument in Earth coordinates.',
        'units': 'uT'
    },
    'magnetometer_y': {
        'long_name': 'Magnetometer Y',
        'comment': 'Magnetic field strength (Hy) in the x direction. Units are not defined in the manual, but are'
                   'believed to be uT. Combined with the pitch and roll sensor data, the magnetometer data is used'
                   'to determine the heading of the instrument in Earth coordinates.',
        'units': 'uT'
    },
    'magnetometer_z': {
        'long_name': 'Magnetometer Z',
        'comment': 'Magnetic field strength (Hz) in the x direction. Units are not defined in the manual, but are'
                   'believed to be uT. Combined with the pitch and roll sensor data, the magnetometer data is used'
                   'to determine the heading of the instrument in Earth coordinates.',
        'units': 'uT'
    },
    'acceleration_x': {
        'long_name': 'Acceleration X',
        'comment': 'Acceleration in the x direction converted from counts to gravity by dividing by 16384. Can be '
                   'further converted to m/s^2 by multiplying by the gravitational constant of 9.80665 m/s^2.',
        'units': 'gravity'
    },
    'acceleration_y': {
        'long_name': 'Acceleration Y',
        'comment': 'Acceleration in the y direction converted from counts to gravity by dividing by 16384. Can be '
                   'further converted to m/s^2 by multiplying by the gravitational constant of 9.80665 m/s^2.',
        'units': 'gravity'
    },
    'acceleration_z': {
        'long_name': 'Acceleration Z',
        'comment': 'Acceleration in the z direction converted from counts to gravity by dividing by 16384. Can be '
                   'further converted to m/s^2 by multiplying by the gravitational constant of 9.80665 m/s^2.',
        'units': 'gravity'
    },
    'ambiguity_velocity': {
        'long_name': 'Ambiguity Velocity',
        'comment': ('Ambiguity Velocity is the maximum beam velocity that can be measured by the instrument in its '
                    'deployed configuration. This is set by defining the velocity range during instrument '
                    'configuration.'),
        'units': 'm s-1'
    },
    'beam_mapping': {
        'long_name': 'Beam Mapping',
        'comment': ('An array of 5 values defining which physical beam (1 thru 4) is assigned to the 1st through the '
                    '3rd data sets (the 5th beam and the 4th and 5th data sets are not used). The Aquadopp II assigns '
                    'physical beams 1, 2 and 4 during an upcast to the 1st, 2nd and 3rd data sets, respectively. '
                    'During a downcast, physical beams 2, 3, and 4 are assigned to the 1st thru 3rd data sets. The '
                    'velocity, amplitude, and correlation parameters are named to correspond to the 1st, 2nd, and 3rd '
                    'data sets. With the mapping provided by this parameter, the user can determine which physical '
                    'beam is assigned to each data set.'),
        # 'units': '',    deliberately left blank, no units for this value
    },
    'beam_map': {
        'long_name': 'Beam Map',
        'comment': 'Indexing array used with the beam_mapping parameter to map the physical beam to the data set.',
        'units': 'count'
    },
    'transmit_energy': {
        'long_name': 'Transmit Energy',
        'comment': 'Transmit energy setting of the Aquadopp II.',
        # 'units': '',    deliberately left blank, don't know the units for this value, could be dB?
    },
    'velocity_scaling': {
        'long_name': 'Velocity Scaling Exponent',
        'comment': ('The velocity scaling exponent, VScale, is used to convert the beam velocity data to m/s. The '
                    'velocity data is scaled by multiplying the reported value by 10^VScale.'),
        # 'units': '',    deliberately left blank, no units for this value
    },
    'power_level': {
        'long_name': 'Power Level',
        'comment': 'Configured power level setting of the Aquadopp II, reported in dB.',
        'units': 'dB'
    },
    'velocity_1': {
        'long_name': 'Velocity Data Set 1',
        'comment': ('The first velocity data set, mapped to a physical beam based on the data set description. Values'
                    'are beam velocities converted to m/s by multiplying the reported value by 10^VScale (the '
                    'velocity scaling exponent).'),
        'data_product_identifier': 'VELPTMN-1ST_L0',
        'units': 'm s-1',
        'ancillary_variables': 'velocity_scaling beam_mapping'
    },
    'velocity_2': {
        'long_name': 'Velocity Data Set 2',
        'comment': ('The second velocity data set, mapped to a physical beam based on the data set description. Values'
                    'are beam velocities converted to m/s by multiplying the reported value by 10^VScale (the '
                    'velocity scaling exponent).'),
        'data_product_identifier': 'VELPTMN-2ND_L0',
        'units': 'm s-1',
        'ancillary_variables': 'velocity_scaling beam_mapping'
    },
    'velocity_3': {
        'long_name': 'Velocity Data Set 3',
        'comment': ('The third velocity data set, mapped to a physical beam based on the data set description. Values'
                    'are beam velocities converted to m/s by multiplying the reported value by 10^VScale (the '
                    'velocity scaling exponent).'),
        'data_product_identifier': 'VELPTMN-3RD_L0',
        'units': 'm s-1',
        'ancillary_variables': 'velocity_scaling beam_mapping'
    },
    'amplitude_1': {
        'long_name': 'Amplitude Data Set 1',
        'comment': ('Raw measurement, for the first data set, of the difference in frequency between the transmitted '
                    'and the received pulse, which is proportional to the velocity of the water.'),
        'units': 'count'
    },
    'amplitude_2': {
        'long_name': 'Amplitude Data Set 2',
        'comment': ('Raw measurement, for the second data set, of the difference in frequency between the transmitted '
                    'and the received pulse, which is proportional to the velocity of the water.'),
        'units': 'count'
    },
    'amplitude_3': {
        'long_name': 'Amplitude Data Set 3',
        'comment': ('Raw measurement, for the third data set, of the difference in frequency between the transmitted '
                    'and the received pulse, which is proportional to the velocity of the water.'),
        'units': 'count'
    },
    'correlation_1': {
        'long_name': 'Percent Correlation Data Set 1',
        'comment': ('Percent correlation is a measure of the similarity of two pulse echoes being measured by the '
                    'Doppler instrument. Zero correlation means nothing at all is similar between the two echoes, '
                    'whereas a correlation of 100 means the two echoes are identical. We want high correlation '
                    'because it gives us confidence the system measured the two pulses it originally sent out and '
                    'is determining a valid phase shift. Values less than 50% should be discarded.'),
        'units': 'percent'
    },
    'correlation_2': {
        'long_name': 'Percent Correlation Data Set 2',
        'comment': ('Percent correlation is a measure of the similarity of two pulse echoes being measured by the '
                    'Doppler instrument. Zero correlation means nothing at all is similar between the two echoes, '
                    'whereas a correlation of 100 means the two echoes are identical. We want high correlation '
                    'because it gives us confidence the system measured the two pulses it originally sent out and '
                    'is determining a valid phase shift. Values less than 50% should be discarded.'),
        'units': 'percent'
    },
    'correlation_3': {
        'long_name': 'Percent Correlation Data Set 3',
        'comment': ('Percent correlation is a measure of the similarity of two pulse echoes being measured by the '
                    'Doppler instrument. Zero correlation means nothing at all is similar between the two echoes, '
                    'whereas a correlation of 100 means the two echoes are identical. We want high correlation '
                    'because it gives us confidence the system measured the two pulses it originally sent out and '
                    'is determining a valid phase shift. Values less than 50% should be discarded.'),
        'units': 'percent'
    },
    'relative_velocity_east': {
        'long_name': 'Eastward Sea Water Velocity',
        'standard_name': 'eastward_sea_water_velocity',
        'comment': ('Eastward sea water velocity corrected for magnetic declination and scaled to m/s. Note, '
                    'this is a relative velocity, not an absolute velocity, as the movement of the profiler '
                    'during ascent/descent has not been removed.'),
        'data_product_identifier': 'VELPTMN-VLE_L1',
        'units': 'm s-1',
        'ancillary_variables': 'beam_mapping velocity_1 velocity_2nd velocity_3rd heading pitch roll'
    },
    'relative_velocity_north': {
        'long_name': 'Northward Sea Water Velocity',
        'standard_name': 'northward_sea_water_velocity',
        'comment': ('Northward sea water velocity corrected for magnetic declination and scaled to m/s. Note, '
                    'this is a relative velocity, not an absolute velocity, as the movement of the profiler '
                    'during ascent/descent has not been removed.'),
        'data_product_identifier': 'VELPTMN-VLN_L1',
        'units': 'm s-1',
        'ancillary_variables': 'beam_mapping velocity_1 velocity_2nd velocity_3rd heading pitch roll'
    },
    'relative_velocity_vertical': {
        'long_name': 'Upward Sea Water Velocity',
        'standard_name': 'upward_sea_water_velocity',
        'comment': ('Vertical sea water velocity component scaled to m/s. Note, this is a relative velocity, not an '
                    'absolute velocity, as the movement of the profiler during ascent/descent has not been removed.'),
        'data_product_identifier': 'VELPTMN-VLU_L1',
        'units': 'm s-1',
        'ancillary_variables': 'beam_mapping velocity_1 velocity_2nd velocity_3rd heading pitch roll'
    },
    'aquadopp_sensor_quality_flag': {
        'long_name': 'Aquadopp II Sensor Quality Flag',
        'standard_name': 'aggregate_quality_flag',
        'comment': ('The Aquadopp II sensor quality flag is built from a set of vendor defined quality checks looking '
                    'at a combination of different variables to create a single automated quality assessment of the '
                    'data. The final flag value represents the worst case assessment of the data quality. Note, this '
                    'flag represents only an automated quality assessment of the data and should be used in '
                    'conjunction with other flags and data quality information.'),
        'flag_values': np.array([1, 2, 3, 4, 9]),
        'flag_meanings': 'pass not_evaluated suspect_or_of_high_interest fail missing'
    }
})


def quality_checks(ds):
    """
    Quality assessment of the pitch, roll, and pressure values for the VEL3D
    using a subset of the QARTOD flags to indicate the quality. QARTOD
    flags used are:

        1 = Pass
        3 = Of High Interest or Suspect
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
    if 'speed_of_sound' in ds.variables:
        m = (ds['speed_of_sound'] < 1400) | (ds['speed_of_sound'] > 1700)
        qc_flag[m] = 4

    # test for pressure out of range (catch those periods when the instrument is being deployed/recovered or
    # out of the water)
    if 'sea_water_pressure' in ds.variables:
        m = ds['sea_water_pressure'] <= 15
        qc_flag[m] = 4

    # test for any error codes, which indicate a problem with the instrument
    if 'error_code' in ds.variables:
        m = (ds['error_code'].astype(int) & 1) == 1
        qc_flag[m] = 4

    # vector: test for low correlation values (less than 50%) from any of the three beams
    if 'correlation_beam1' in ds.variables:
        m = (ds['correlation_beam1'] < 50) | (ds['correlation_beam2'] < 50) | (ds['correlation_beam3'] < 50)
        qc_flag[m] = 4

    # aquadopp: test for low correlation values (less than 50%) from any of the three data sets
    if 'correlation_1' in ds.variables:
        m = (ds['correlation_1'] < 50) | (ds['correlation_2'] < 50) | (ds['correlation_3'] < 50)
        qc_flag[m] = 4

    # aquadopp: test for velocity values greater than the ambiguity velocity
    if 'ambiguity_velocity' in ds.variables:
        m = (np.abs(ds['velocity_1']) > ds['ambiguity_velocity']) | \
            (np.abs(ds['velocity_2']) > ds['ambiguity_velocity']) | \
            (np.abs(ds['velocity_3']) > ds['ambiguity_velocity'])
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
        a burst of data. Used for QC testing.
    :param system: xarray dataset with the system data. This data set contains
        the 1 Hz instrument status data (e.g. battery voltage, temperature,
        pressure, tilt, roll, etc.). Used for QC testing and for environmental
        context.
    :param velocity: xarray dataset with the velocity data. This data set
        contains the 8 Hz velocity data collected by the instrument.
    :param burst: boolean, used to indicate if the data should be burst
        averaged or not. Default is False.
    :return ds: cleaned up data set
    """
    # drop some of the variables: velocity
    #   internal_timestamp == time == redundant, so can remove
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
    #   turbulent_velocity_vertical == vel3d_c_upward_turbulent_velocity == redundant but properly scaled
    #   vel3d_c_upward_turbulent_velocity_qc_executed == drop in favor of instrument specific tests
    #   vel3d_c_upward_turbulent_velocity_qc_results == drop in favor of instrument specific tests
    drop_vars = ['internal_timestamp', 'analog_input_1', 'analog_input_2', 'sea_water_pressure_mbar_qc_executed',
                 'sea_water_pressure_mbar_qc_results', 'turbulent_velocity_east_qc_executed',
                 'turbulent_velocity_east_qc_results', 'turbulent_velocity_north_qc_executed',
                 'turbulent_velocity_north_qc_results', 'turbulent_velocity_vertical_qc_executed',
                 'turbulent_velocity_vertical_qc_results', 'vel3d_c_eastward_turbulent_velocity_qc_executed',
                 'vel3d_c_eastward_turbulent_velocity_qc_results', 'vel3d_c_northward_turbulent_velocity_qc_executed',
                 'vel3d_c_northward_turbulent_velocity_qc_results', 'turbulent_velocity_vertical',
                 'vel3d_c_upward_turbulent_velocity_qc_executed', 'vel3d_c_upward_turbulent_velocity_qc_results'
                 ]
    for var in velocity.variables:
        if var in drop_vars:
            velocity = velocity.drop_vars(var)

    # drop some of the variables: system
    #   internal_timestamp == time == redundant, so can remove
    #   analog_input == not used, no data
    #   date_time_string == time == redundant, so can remove
    #   speed_of_sound_qc_executed == drop in favor of instrument specific tests
    #   speed_of_sound_qc_results == drop in favor of instrument specific tests
    drop_vars = ['internal_timestamp', 'analog_input', 'date_time_string',
                 'speed_of_sound_qc_executed', 'speed_of_sound_qc_results']
    for var in system.variables:
        if var in drop_vars:
            system = system.drop_vars(var)

    # drop some of the variables: header
    #   internal_timestamp == time == redundant, so can remove
    #   date_time_string == redundant, so can remove
    #   number_velocity_records == this is always 1440, so not useful
    drop_vars = ['internal_timestamp', 'date_time_string', 'number_velocity_records']
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
        'temperature_centidegree': 'sea_water_temperature',
        # velocity packets
        'sea_water_pressure_mbar': 'sea_water_pressure',  # mbar units are incorrect, should be 0.001 dbar
        'amplitude_beam_1': 'amplitude_beam1',
        'amplitude_beam_2': 'amplitude_beam2',
        'amplitude_beam_3': 'amplitude_beam3',
        'correlation_beam_1': 'correlation_beam1',
        'correlation_beam_2': 'correlation_beam2',
        'correlation_beam_3': 'correlation_beam3',
        'turbulent_velocity_east': 'velocity_east',
        'turbulent_velocity_north': 'velocity_north',
        'vel3d_c_eastward_turbulent_velocity': 'velocity_vertical',
        'vel3d_c_eastward_turbulent_velocity': 'velocity_east_corrected',
        'vel3d_c_northward_turbulent_velocity': 'velocity_north_corrected',
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
    system['sea_water_temperature'] = system['sea_water_temperature'] / 100.0  # convert from cdegC to degC
    system['speed_of_sound'] = system['speed_of_sound'] / 10.0  # convert from dm/s to m/s
    velocity['sea_water_pressure'] = velocity['sea_water_pressure'] / 1000.0  # convert from 0.001 dbar to dbar

    # pull the scaling factor out of the status code
    status_code = system.status_code.values
    scaling = np.array([[n >> i & 1 for i in range(0, int(n).bit_length())] for n in status_code])[:, 1]
    system['scaling_factor'] = ('time', np.where(scaling == 1, 0.1, 1))

    # merge the three data sets together and then use a forward-fill to match header and system data
    # with the velocity data
    velocity['time'] = velocity['time'] - np.timedelta64(2, 's')  # adjust incorrect time offset of 2 seconds
    system['time'] = system['time'] - np.timedelta64(1, 's')  # adjust incorrect time offset of 1 second
    vel3d = xr.merge([header, system, velocity])
    vel3d = vel3d.ffill(dim='time')

    # adjust the uncorrected velocity data based on the scaling factor and convert to m/s
    vel3d['velocity_east'] = vel3d.velocity_east * vel3d.scaling_factor / 1000.0
    vel3d['velocity_north'] = vel3d.velocity_north * vel3d.scaling_factor / 1000.0

    # add the SNR values to the data set (useful for QC purposes, though correlation percentages are a better metric)
    vel3d['snr_beam1'] = (vel3d['amplitude_beam1'] - vel3d['noise_amplitude_beam1']) * 0.43
    vel3d['snr_beam2'] = (vel3d['amplitude_beam2'] - vel3d['noise_amplitude_beam2']) * 0.43
    vel3d['snr_beam3'] = (vel3d['amplitude_beam3'] - vel3d['noise_amplitude_beam3']) * 0.43

    # now add QC flags to the data set
    vel3d['vector_sensor_quality_flag'] = quality_checks(vel3d)

    # reset some attributes
    for key, value in VECTOR.items():
        for atk, atv in value.items():
            if key in vel3d.variables:
                vel3d[key].attrs[atk] = atv

    # add the original variable name as an attribute, if renamed
    for key, value in rename.items():
        vel3d[value].attrs['ooinet_variable_name'] = key

    # burst averaging the data, if requested
    if burst:
        # use the quality flag to remove bad data prior to burst averaging
        m = vel3d['vector_sensor_quality_flag'] == 3
        vel3d['velocity_vertical'] = vel3d['velocity_vertical'].where(~m)
        vel3d['velocity_east'] = vel3d['velocity_east'].where(~m)
        vel3d['velocity_east_corrected'] = vel3d['velocity_east_corrected'].where(~m)
        vel3d['velocity_north'] = vel3d['velocity_north'].where(~m)
        vel3d['velocity_north_corrected'] = vel3d['velocity_north_corrected'].where(~m)

        # convert the heading to radians before burst averaging
        vel3d['heading'] = np.unwrap(np.deg2rad(vel3d['heading']))

        # resample the data to 30 minute intervals using a median
        vel3d = vel3d.resample(time='1800s').median(dim='time', keep_attrs=True, skipna=True)
        vel3d = vel3d.where(~np.isnan(vel3d.deployment), drop=True)  # drop the fill values

        # convert the heading back to degrees after burst averaging
        vel3d['heading'] = np.mod(np.rad2deg(vel3d['heading']), 360)

    return vel3d


def mmp_aquadopp(ds, binning=False, bin_size=2.0):
    """
    Takes Nortek Aquadopp II data, recorded by the McLane Moored Profiler (MMP)
    and cleans up the data set to make it more user-friendly. Primary task is
    re-working the beam amplitude, velocity and correlation data so it is more
    understandable and user-friendly, and adding missing variables to the data
    depending on whether this is a decimated data set (telemetered) or not
    (recovered) so full deployments can be cross-compared.

    :param ds: xarray dataset with the MMP Aquadopp data
    :return ds: cleaned up data set, fully populated with all the necessary
        variables for further analysis of the Aquadopp II data.
    """
    # drop some of the variables
    #   internal_timestamp == time == redundant, so can remove
    #   vel3d_k_version == always set to the same value, so not useful
    #   vel3d_k_id == Nortek binary data packet ID, so not useful
    #   vel3d_k_beams == always set to the same value, so not useful
    #   vel3d_k_configuration == always set to the same configuration, so not useful
    #   vel3d_k_beams_coordinate == always set to the same value, so not useful
    #   vel3d_k_serial == serial number is stored in the global attributes
    #   vel3d_k_pressure_qc_executed == test limits did NOT account for the scaling of 0.001 dbar
    #   vel3d_k_pressure_qc_results == test limits did NOT account for the scaling of 0.001 dbar
    drop_vars = ['internal_timestamp', 'vel3d_k_version', 'vel3d_k_id', 'vel3d_k_beams', 'vel3d_k_configuration',
                 'vel3d_k_beams_coordinate', 'vel3d_k_serial', 'vel3d_k_pressure_qc_executed',
                 'vel3d_k_pressure_qc_results']
    for var in ds.variables:
        if var in drop_vars:
            ds = ds.drop_vars(var)

    # The original calculation of the time failed to include the microsecond values in the date_time_array for
    # the recovered data. Indeed, the original calculation assumes a fixed data rate rather than just using the
    # already existing microseconds values because......not sure why. Re-calculate the time using the date_time_array
    # and the vel3d_k_micro_second variable and then remove the date_time_array variable and the dimension variable.
    if 'vel3d_k_micro_second' in ds.variables:
        usec = ds['vel3d_k_micro_second'].values.astype(int) * 100  # convert from 100 microseconds to microseconds
        dt_array = ds['date_time_array'].values.astype(int)
        dt_array = np.append(dt_array, np.atleast_2d(usec).T, axis=1)
        dt_array[:, 0] = dt_array[:, 0] + 1900  # convert from years since 1900 to years since 0000
        dt_array[:, 1] = dt_array[:, 1] + 1  # convert from 0-based to 1-based months
        dt = np.zeros(dt_array.shape[0], dtype='datetime64[ns]')
        for i in range(dt_array.shape[0]):
            dt[i] = np.datetime64('%04d-%02d-%02dT%02d:%02d:%02d.%06d' % tuple(dt_array[i]))

        # reset the time variable and remove the vel3d_k_micro_second variable
        ds['time'] = dt
        ds = ds.drop_vars(['vel3d_k_micro_second'])

    # drop the date_time_array variable and the dimension variable
    ds = ds.drop_vars(['date_time_array', 'date_time_array_dim_0'])

    # make sure our time record is monotonically increasing
    _, index = np.unique(ds['time'], return_index=True)
    ds = ds.isel(time=index)

    # rename some parameters here to get a better defined data set with cleaner attributes
    rename = {
        'vel3d_k_speed_sound': 'speed_of_sound',
        'vel3d_k_temp_c': 'sea_water_temperature',
        'int_ctd_pressure': 'ctd_pressure',
        'vel3d_k_pressure': 'sea_water_pressure',
        'vel3d_k_heading': 'heading',
        'vel3d_k_pitch': 'pitch',
        'vel3d_k_roll': 'roll',
        'vel3d_k_transmit_energy': 'transmit_energy',
        'vel3d_k_error': 'error_code',
        'vel3d_k_status': 'status_code',
        'vel3d_k_cell_size': 'cell_size',
        'vel3d_k_blanking': 'blanking_distance',
        'vel3d_k_velocity_range': 'velocity_range',
        'vel3d_k_battery_voltage': 'battery_voltage',
        'vel3d_k_mag_x': 'magnetometer_x',
        'vel3d_k_mag_y': 'magnetometer_y',
        'vel3d_k_mag_z': 'mahnetometer_z',
        'vel3d_k_acc_x': 'acceleration_x',
        'vel3d_k_acc_y': 'acceleration_y',
        'vel3d_k_acc_z': 'acceleration_z',
        'vel3d_k_ambiguity': 'ambiguity_velocity',
        'vel3d_k_data_set_description': 'beam_mapping',
        'vel3d_k_data_set_description_dim_0': 'beam_map',
        'vel3d_k_transmit_energy': 'transmit_energy',
        'vel3d_k_v_scale': 'velocity_scaling',
        'vel3d_k_power_level': 'power_level',
        'vel3d_k_vel0': 'velocity_1',
        'vel3d_k_vel1': 'velocity_2',
        'vel3d_k_vel2': 'velocity_3',
        'vel3d_k_amp1': 'amplitude_1',
        'vel3d_k_amp2': 'amplitude_2',
        'vel3d_k_amp0': 'amplitude_3',
        'vel3d_k_corr0': 'correlation_1',
        'vel3d_k_corr1': 'correlation_2',
        'vel3d_k_corr2': 'correlation_3',
        'vel3d_k_eastward_velocity': 'relative_velocity_east',
        'vel3d_k_eastward_velocity_qc_executed': 'relative_velocity_east_qc_executed',
        'vel3d_k_eastward_velocity_qc_results': 'relative_velocity_east_qc_results',
        'vel3d_k_northward_velocity': 'relative_velocity_north',
        'vel3d_k_northward_velocity_qc_executed': 'relative_velocity_north_qc_executed',
        'vel3d_k_northward_velocity_qc_results': 'relative_velocity_north_qc_results',
        'vel3d_k_upward_velocity': 'relative_velocity_vertical',
        'vel3d_k_upward_velocity_qc_executed': 'relative_velocity_vertical_qc_executed',
        'vel3d_k_upward_velocity_qc_results': 'relative_velocity_vertical_qc_results',
    }
    for key in rename.keys():
        if key in ds.variables:
            ds = ds.rename({key: rename.get(key)})

    # convert some variables to more standard units
    ds['heading'] = ds['heading'] * 0.1  # convert from ddeg to deg
    ds['pitch'] = ds['pitch'] * 0.1  # convert from ddeg to deg
    ds['roll'] = ds['roll'] * 0.1  # convert from ddeg to deg
    ds['sea_water_temperature'] = ds['sea_water_temperature'] * 0.01  # convert from cdegC to degC
    ds['velocity_1'] = ds['velocity_1'] * 10. ** ds['velocity_scaling']  # convert to m/s
    ds['velocity_2'] = ds['velocity_2'] * 10. ** ds['velocity_scaling']  # convert to m/s
    ds['velocity_3'] = ds['velocity_3'] * 10. ** ds['velocity_scaling']  # convert to m/s

    if 'vel3d_k_wfp_instrument' in ds.attrs['source']:  # recovered data has the following additional variables
        ds['sea_water_pressure'] = ds['sea_water_pressure'] * 0.001  # convert from 0.001 dbar to dbar
        ds['battery_voltage'] = ds['battery_voltage'] * 0.1  # convert from dV to V
        ds['speed_of_sound'] = ds['speed_of_sound'] * 0.1  # convert from dm/s to m/s
        ds['ambiguity_velocity'] = ds['ambiguity_velocity'] * 0.1 / 1000 # convert from 0.1 mm/s to m/s
        ds['acceleration_x'] = ds['acceleration_x'] / 16384.  # convert from raw counts to gravity
        ds['acceleration_y'] = ds['acceleration_y'] / 16384.  # convert from raw counts to gravity
        ds['acceleration_z'] = ds['acceleration_z'] / 16384.  # convert from raw counts to gravity

    # add a profile id to the data set to enable easy separation of the profiles
    ds = create_profile_id(ds)

    # now add QC flags to the data set
    ds['aquadopp_sensor_quality_flag'] = quality_checks(ds)

    # reset some attributes
    for key, value in AQUADOPP.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # add the original variable name as an attribute, if renamed
    for key, value in rename.items():
        if value in ds.variables:
            ds[value].attrs['ooinet_variable_name'] = key

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, of interest == 3, and fail == 4.
    ds = parse_qc(ds)

    if binning:
        # use the quality flag to remove bad data prior to depth binning
        m = ds['aquadopp_sensor_quality_flag'] == 4
        ds['relative_velocity_east'] = ds['relative_velocity_east'].where(~m)
        ds['relative_velocity_north'] = ds['relative_velocity_north'].where(~m)
        ds['relative_velocity_vertical'] = ds['relative_velocity_vertical'].where(~m)

        # convert the heading to radians before binning
        ds['heading'] = np.deg2rad(ds['heading'])

        # group the data by the profile number and then bin the data into 2 m depth bins
        # (nominal ascent rate of the MMP is 15-25 cm/s)
        vocab = get_vocabulary(ds.attrs['subsite'], ds.attrs['node'], ds.attrs['sensor'])[0]
        site_depth = vocab['maxdepth'] - 20  # ~20 meters from the bottom
        profiles = ds.groupby('profile')
        profiles = [profile[1] for profile in profiles]
        partial_binning = partial(bin_profiles, site_depth=site_depth, bin_size=bin_size)
        with ProcessPoolExecutor(max_workers=N_CORES) as executor:
            binned = list(tqdm(executor.map(partial_binning, profiles), total=len(profiles),
                               desc='Smoothing and binning each profile into 2 m depth bins', file=sys.stdout))

        # reset the dataset now using binned profiles
        binned = [i for i in binned if i is not None]
        binned = xr.concat(binned, 'time')
        binned = binned.sortby(['deployment', 'profile', 'time'])

        # convert the heading back to degrees after binning
        binned['heading'] = np.mod(np.rad2deg(binned['heading']), 360)

        # make sure our time record is monotonically increasing
        _, index = np.unique(binned['time'], return_index=True)
        binned = binned.isel(time=index)

        # reset the original integer variables to integers after the binning
        for v in binned.variables:
            if np.issubdtype(ds[v].dtype, np.integer):
                binned[v] = binned[v].astype(ds[v].dtype)

        # reset the dataset with the new binned profiles
        ds = binned.copy()

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
        tag = '.*deployment%04d.*VEL3D.*\\.nc$' % deploy
        velocity = load_gc_thredds(site, node, sensor, method, stream, tag)

        # get the deployment dates, which will be used to request the system and header data from the M2M system
        start, stop = get_deployment_dates(site, node, sensor, deploy)

        # check to see if we downloaded any data
        if not velocity:
            exit_text = ('Velocity data unavailable for %s-%s-%s, %s, %s, deployment %d.' % (site, node, sensor, method,
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
        tag = '.*VEL3D.*\\.nc$'
        velocity = m2m_collect(r, tag)

        # check to see if we downloaded any data
        if not velocity:
            exit_text = ('Velocity data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                             stream, start, stop))
            raise SystemExit(exit_text)

    if node != 'WFP01':
        # we are working with the Nortek Vector, so we need to get the system and header data from the M2M system.
        # the additional system and header data are only available from the M2M system because they are not
        # considered "science" data, but rather "engineering" data which is not included in the GC THREDDS catalog.
        system_stream = re.sub('velocity', 'system', stream)
        r = m2m_request(site, node, sensor, method, system_stream, start, stop)
        if not r:
            exit_text = ('Request failed for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                  system_stream, start, stop))
            raise SystemExit(exit_text)
        system = m2m_collect(r, tag)  # provides heading, pitch, roll, and seawater temperature data
        if not system:
            exit_text = ('System data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                           system_stream, start, stop))
            raise SystemExit(exit_text)

        header_stream = re.sub('velocity', 'header', stream)
        r = m2m_request(site, node, sensor, method, header_stream, start, stop)
        if not r:
            exit_text = ('Request failed for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                  header_stream, start, stop))
            raise SystemExit(exit_text)
        header = m2m_collect(r, tag)  # provides noise floor measurements needed to QC the data
        if not system:
            exit_text = ('System data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                           header_stream, start, stop))
            raise SystemExit(exit_text)

        # clean-up and reorganize the data
        vel3d = vel3d_datalogger(header, system, velocity, burst=burst)
        depth = vel3d.depth.mean().values
        vel3d = update_dataset(vel3d, depth)
    else:
        # we are working with the Nortek Aquadopp II on the McLane Moored Profiler (MMP), so we just need to clean up
        # and reorganize the data to make it more user-friendly
        vel3d = mmp_aquadopp(velocity, binning=burst, bin_size=2.0)

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    vel3d.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
