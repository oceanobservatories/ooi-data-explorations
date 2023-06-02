#!/usr/bin/env python
# -*- coding: utf-8 -*-
import dask
import numpy as np
import os
import re
import sys
import xarray as xr

from copy import copy
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from dask.diagnostics import ProgressBar
from scipy.interpolate import CubicSpline
from tqdm import tqdm

from ooi_data_explorations.calibrations import Coefficients
from ooi_data_explorations.common import inputs, get_vocabulary, get_calibrations_by_uid, m2m_request, \
    list_files, m2m_collect, load_gc_thredds, update_dataset, ENCODINGS
from ooi_data_explorations.profilers import create_profile_id, bin_profiles

from pyseas.data.opt_functions import opt_internal_temp, opt_external_temp
from pyseas.data.opt_functions_tscor import tscor

# reset the variable level attributes and set some global defaults
N_CORES = int(os.cpu_count() / 2) - 1
N_WORKERS = int(os.cpu_count() / 4) - 1
FILL_INT = -9999999
ATTRS = dict({
    # parsed (raw) variables and attributes
    'serial_number': {
        'long_name': 'Unit Serial Number',
        # 'units': ''    # deliberately left blank, no units for this value
        'comment': ('Instrument serial number, useful for tracking history of instrument through deployments and '
                    'to make sure the correct calibration coefficients have been applied.')
    },
    'external_temp_raw': {
        'long_name': 'Raw In-Situ Temperature',
        'units': 'count',
        'data_product_identifier': 'OPTTEMP_L0',
        'comment': ('Raw measurements, reported in counts, from the AC-S external temperature sensor. This sensor '
                    'measures the in-situ seawater temperature.')
    },
    'internal_temp_raw': {
        'long_name': 'Raw Internal Instrument Temperature',
        'units': 'count',
        'comment': ('Raw measurements, reported in counts, from the AC-S internal temperature sensor. This sensor '
                    'measures the internal instrument temperature and is used in converting the raw optical '
                    'measurements into absorbance and attenuation estimates.')
    },
    'elapsed_run_time': {
        'long_name': 'Elapsed Run Time',
        'units': 'ms',
        'comment': 'Time in milliseconds since the instrument was powered on.'
    },
    'wavelength_number': {
        'long_name': 'Wavelength Number',
        'units': 'count',
        'comment': ('An index counter between 0 and 99 used to set a common length dimension for the absorbance and '
                    'attenuation measurements. The actual number of wavelengths is variable between sensors '
                    'and may even change for a particular sensor over time if servicing requires a replacement '
                    'of the filter set. The actual number of wavelengths for this sensor is represented here '
                    'by the attribute actual_wavelengths.')
        # 'actual_wavelengths': ''  # deliberately left blank, created during the processing
    },
    'wavelength_a': {
        'long_name': 'Absorption Channel Wavelengths',
        'standard_name': 'radiation_wavelength',
        'units': 'nm',
        'comment': ('Absorbance channel measurement wavelengths, specific to the filter wheel set installed in '
                    'the AC-S.'),
        '_FillValue': np.nan
    },
    'a_reference_dark': {
        'long_name': 'Absorption Channel Dark Reference',
        'units': 'count',
        'comment': ('Absorption channel reference detector dark counts (before the lamp is turned on). Used in '
                    'conversion of the raw absorption channel measurements to absorbance estimates.')
    },
    'a_reference': {
        'long_name': 'Absorption Channel Reference Measurements',
        'units': 'count',
        'comment': ('Absorption channel reference detector raw counts (while the lamp is turned on). Used in '
                    'conversion of the raw absorption channel measurements to absorbance estimates.'),
        'data_product_identifier': 'OPTAREF_L0',
        '_FillValue': FILL_INT
    },
    'a_signal_dark': {
        'long_name': 'Absorption Channel Dark Signal',
        'units': 'count',
        'comment': ('Absorption channel signal detector dark counts (before the lamp is turned on). Used in conversion '
                    'of the raw absorption channel measurements to absorbance estimates.')
    },
    'a_signal': {
        'long_name': 'Absorption Channel Signal Measurements',
        'units': 'count',
        'comment': ('Absorption channel signal detector raw counts (while the lamp is turned on). Used in conversion '
                    'of the raw absorption channel measurements to absorbance estimates.'),
        'data_product_identifier': 'OPTASIG_L0',
        '_FillValue': FILL_INT
    },
    'wavelength_c': {
        'long_name': 'Attenuation Channel Wavelengths',
        'standard_name': 'radiation_wavelength',
        'units': 'nm',
        'comment': ('Attenuation channel measurement wavelengths, specific to the filter wheel set installed in '
                    'the AC-S.'),
        '_FillValue': np.nan
    },
    'c_reference_dark': {
        'long_name': 'Attenuation Channel Dark Reference',
        'units': 'count',
        'comment': ('Attenuation channel reference detector dark counts (before the lamp is turned on). Used in '
                    'conversion of the raw attenuation channel measurements to attenuation estimates.')
    },
    'c_reference': {
        'long_name': 'Attenuation Channel Reference Measurements',
        'units': 'count',
        'comment': ('Attenuation channel reference detector raw counts (while the lamp is turned on). Used in '
                    'conversion of the raw attenuation channel measurements to attenuation estimates.'),
        'data_product_identifier': 'OPTCREF_L0',
        '_FillValue': FILL_INT
    },
    'c_signal_dark': {
        'long_name': 'Attenuation Channel Dark Signal',
        'units': 'count',
        'comment': ('Attenuation channel signal detector dark counts (before the lamp is turned on). Used in '
                    'conversion of the raw attenuation channel measurements to attenuation estimates.')
    },
    'c_signal': {
        'long_name': 'Attenuation Channel Signal Measurements',
        'units': 'count',
        'comment': ('Attenuation channel signal detector raw counts (while the lamp is turned on). Used in conversion '
                    'of the raw attenuation channel measurements to attenuation estimates.'),
        'data_product_identifier': 'OPTCSIG_L0',
        '_FillValue': FILL_INT
    },

    # Data from a co-located CTD, if available, interpolated into the data set
    'sea_water_temperature': {
        'long_name': 'Seawater Temperature',
        'standard_name': 'sea_water_temperature',
        'units': 'degree_Celsius',
        'comment': ('Sea water temperature is the in situ temperature of the sea water near the sensor. Measurements '
                    'are from a co-located CTD'),
        'data_product_identifier': 'TEMPWAT_L1',
        '_FillValue': np.nan
    },
    'sea_water_practical_salinity': {
        'long_name': 'Practical Salinity',
        'standard_name': 'sea_water_practical_salinity',
        'units': '1',
        'comment': ('Salinity is generally defined as the concentration of dissolved salt in a parcel of sea water. '
                    'Practical Salinity is a more specific unit-less quantity calculated from the conductivity of '
                    'sea water and adjusted for temperature and pressure. It is approximately equivalent to Absolute '
                    'Salinity (the mass fraction of dissolved salt in sea water), but they are not interchangeable. '
                    'Measurements are from a co-located CTD.'),
        'data_product_identifier': 'PRACSAL_L2',
        '_FillValue': np.nan
    },
    'sea_water_pressure': {
        'long_name': 'Seawater Pressure',
        'standard_name': 'sea_water_pressure',
        'units': 'dbar',
        'comment': ('Seawater Pressure refers to the pressure exerted on a sensor in situ by the weight of the column '
                    'of seawater above it. It is calculated by subtracting one standard atmosphere from the absolute '
                    'pressure at the sensor to remove the weight of the atmosphere on top of the water column. The '
                    'pressure at a sensor in situ provides a metric of the depth of that sensor. Measurements are from '
                    'a co-located CTD'),
        'data_product_identifier': 'PRESWAT_L1',
        '_FillValue': np.nan
    },

    # OOI calculated beam attenuation and optical absorption
    'beam_attenuation': {
        'long_name': 'Optical Beam Attenuation Coefficient',
        'units': 'm-1',
        'comment': ('The Optical Beam Attenuation Coefficient is the rate that the intensity of a beam of light will '
                    'decrease in response to the combined effects of absorption and scatter as a function of '
                    'propagation distance. The Attenuation Coefficient results from the spectral beam attenuation of '
                    'the combination of all seawater impurities including all particulate and dissolved matter of '
                    'optical importance.'),
        'data_product_identifier': 'OPTATTN_L2',
        'ancillary_variables': ('wavelength_c internal_temp c_signal c_reference c_signal_dark c_reference_dark '
                                'seawater_temperature practical_salinity'),
        '_FillValue': np.nan
    },
    'optical_absorption': {
        'long_name': 'Optical Absorption Coefficient',
        'standard_name': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water',
        'units': 'm-1',
        'comment': ('Optical Absorption Coefficient is the rate that the intensity of a beam of light will decrease '
                    'in response to the absorption (removal) of light energy as a function of propagation distance. '
                    'The Optical Absorption Coefficient reflects the absorption coefficient for the combination of all '
                    'seawater impurities including all particulate and dissolved matter of optical importance.'),
        'data_product_identifier': 'OPTABSN_L2',
        'ancillary_variables': ('wavelength_a internal_temp a_signal a_reference a_signal_dark a_reference_dark '
                                'seawater_temperature practical_salinity'),
        '_FillValue': np.nan
    },

    # Derived values in the re-processed data set
    'external_temp': {
        'long_name': 'External Instrument Temperature',
        'standard_name': 'sea_water_temperature',
        'units': 'degrees_Celsius',
        'comment': ('In-situ sea water temperature measurements from the sensor mounted at the top of the '
                    'AC-S pressure housing.'),
        'ancillary_variables': 'external_temp_raw'
    },
    'internal_temp': {
        'long_name': 'Internal Instrument Temperature',
        'units': 'degrees_Celsius',
        'comment': 'Internal instrument temperature, used to convert raw absorbance and attenuation measurements.',
        'ancillary_variables': 'internal_temp_raw'
    },
    'apg': {
        'long_name': 'Particulate and Dissolved Absorbance',
        'units': 'm-1',
        'comment': ('The optical absorption coefficient is the rate that the intensity of a beam of light will '
                    'decrease in response to the absorption (removal) of light energy as a function of propagation '
                    'distance. The optical absorption coefficient reflects the absorption coefficient for the '
                    'combination of all seawater impurities, including all particulate and dissolved matter of '
                    'optical importance.'),
        'data_product_identifier': 'OPTABSN_L1',
        'ancillary_variables': ('wavelength_a internal_temp a_signal_raw a_reference_raw '
                                'a_signal_dark a_reference_dark'),
        '_FillValue': np.nan
    },
    'a_jump_offsets': {
        'long_name': 'Absorption Channel Holographic Grater Jump Offset',
        'units': 'm-1',
        'comment': ('Offset used to correct for spectral jumps commonly seen in the AC-S data where the sensor uses '
                    'two holographic gratings to span the full spectral range. Adding the offset to all values '
                    'from the grate_index (included as an additional attribute) to the end of the spectra will '
                    'restore the AC-S data to values reported by the sensor.'),
        'grate_index': 0,  # Will update in the processing script
        'ancillary_variables': 'wavelength_a apg',
    },
    'apg_ts': {
        'long_name': 'Particulate and Dissolved Absorbance with TS Correction',
        'units': 'm-1',
        'comment': ('The optical absorption coefficient corrected for the effects of temperature and salinity. '
                    'Utilizes data from a co-located CTD for the temperature and salinity, if available. If no '
                    'co-located CTD data is available, will assume a constant salinity of 33 psu and will use '
                    'the OPTAA''s external temperature sensor.'),
        'ancillary_variables': 'wavelength_a sea_water_temperature sea_water_practical_salinity external_temp apg',
        '_FillValue': np.nan
    },
    'apg_ts_s': {
        'long_name': 'Particulate and Dissolved Absorbance with TS and Scatter Correction',
        'standard_name': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water',
        'units': 'm-1',
        'comment': ('The optical absorption coefficient corrected for scattering after correcting for temperature and '
                    'salinity. Utilizes Method 1 (baseline correction) by subtracting the absorption at 715 nm from '
                    'all values; assumes scattering is flat across wavelengths.'),
        'data_product_identifier': 'OPTABSN_L2',
        'ancillary_variables': 'wavelength_a apg_ts',
        '_FillValue': np.nan
    },
    'cpg': {
        'long_name': 'Particulate and Dissolved Attenuation',
        'units': 'm-1',
        'comment': ('The optical beam attenuation coefficient is the rate that the intensity of a beam of light will '
                    'decrease in response to the combined effects of absorption and scatter as a function of '
                    'propagation distance. The attenuation coefficient results from the spectral beam attenuation of '
                    'the combination of all seawater impurities including all particulate and dissolved matter of '
                    'optical importance.'),
        'data_product_identifier': 'OPTATTN_L1',
        'ancillary_variables': ('wavelength_c internal_temp c_signal_raw c_reference_raw '
                                'c_signal_dark c_reference_dark'),
        '_FillValue': np.nan
    },
    'c_jump_offsets': {
        'long_name': 'Attenuation Channel Filter Offsets',
        'units': 'm-1',
        'comment': ('Offset used to correct for spectral jumps commonly seen in the AC-S data where the sensor uses '
                    'two holographic gratings to span the full spectral range. Adding the offset to all values '
                    'from the grate_index (included as an additional attribute) to the end of the spectra will '
                    'restore the AC-S data to values reported by the sensor.'),
        'grate_index': 0,  # Will update in the processing script
        'ancillary_variables': 'wavelength_c cpg',
    },
    'cpg_ts': {
        'long_name': 'Particulate and Dissolved Attenuation with TS Correction',
        'units': 'm-1',
        'comment': ('The optical beam attenuation coefficient corrected for the effects of temperature and salinity. '
                    'Utilizes data from a co-located CTD for the temperature and salinity, if available. If no '
                    'co-located CTD data is available, will assume a constant salinity of 33 psu and will use '
                    'the OPTAA''s external temperature sensor.'),
        'data_product_identifier': 'OPTATTN_L2',
        'ancillary_variables': 'wavelength_c sea_water_temperature sea_water_practical_salinity external_temp cpg',
        '_FillValue': np.nan
    },
    'estimated_chlorophyll': {
        'long_name': 'Estimated Chlorophyll Concentration',
        'standard_name': 'mass_concentration_of_chlorophyll_in_sea_water',
        'units': 'ug L-1',
        'comment': ('Uses the absorption line height at 676 nm, above a linear background between 650 and 715 nm, with '
                    'a chlorophyll specific absorption of 0.020 L/ug/m to estimate the concentration of chlorophyll. '
                    'This method has been shown to be significantly related to extracted chlorophyll concentrations '
                    'and is robust in response to mild to moderate bio-fouling.'),
        'ancillary_variables': 'apg_ts_s',
        '_FillValue': np.nan
    },
    'estimated_poc': {
        'long_name': 'Estimated POC Concentration',
        'standard_name': 'mass_concentration_of_organic_detritus_expressed_as_carbon_in_sea_water',
        'units': 'ug L-1',
        'comment': ('Uses the particulate beam attenuation coefficient at 660 nm and a coefficient of 380 ug/L/m. This '
                    'calculation is not robust in response to bio-fouling and is expected to breakdown as bio-fouling '
                    'begins to dominate the signal.'),
        'ancillary_variables': 'cpg_ts',
        '_FillValue': np.nan
    },
    'ratio_cdom': {
        'long_name': 'CDOM to Chlorophyll Absorbance Ratio',
        'units': '1',
        'comment': ('Ratio of CDOM absorption in the violet portion of the spectrum at 412 nm relative to '
                    'chlorophyll absorption at 440 nm. Ratios greater than 1 indicate a preponderance of CDOM '
                    'absorption relative to chlorophyll.'),
        'ancillary_variables': 'apg_ts_s',
        '_FillValue': np.nan
    },
    'ratio_carotenoids': {
        'long_name': 'Carotenoid to Chlorophyll Absorbance Ratio',
        'units': '1',
        'comment': ('Ratio of carotenoid absorption in the blue-green portion of the spectrum at 490 nm relative to '
                    'chlorophyll absorption at 440 nm. A changing carotenoid to chlorophyll ratio may indicate a shift '
                    'in phytoplankton community composition in addition to changes in light history or bloom health '
                    'and age.'),
        'ancillary_variables': 'apg_ts_s',
        '_FillValue': np.nan
    },
    'ratio_phycobilins': {
        'long_name': 'Phycobilins to Chlorophyll Absorbance Ratio',
        'units': '1',
        'comment': ('Ratio of phycobilin absorption in the green portion of the spectrum at 530 nm relative to '
                    'chlorophyll absorption at 440 nm. Different phytoplankton, notably cyanobacteria, utilize '
                    'phycobilins as accessory light harvesting pigments. An increasing phycobilin to chlorophyll ratio '
                    'may indicate a shift in phytoplankton community composition.'),
        'ancillary_variables': 'apg_ts_s',
        '_FillValue': np.nan
    },
    'ratio_qband': {
        'long_name': 'Chlorophyll Q Band to Soret Band Absorbance Ratio',
        'units': '1',
        'comment': ('The Soret and the Q bands represent the two main absorption bands of chlorophyll. The former '
                    'covers absorption in the blue region of the spectrum, while the latter covers absorption in the '
                    'red region. A decrease in the ratio of the intensity of the Soret band at 440 nm to that of the Q '
                    'band at 676 nm may indicate a change in phytoplankton community structure. All phytoplankton '
                    'contain chlorophyll a as the primary light harvesting pigment, but green algae and '
                    'dinoflagellates contain chlorophyll b and c, respectively, which are spectrally redshifted '
                    'compared to chlorophyll a.'),
        'ancillary_variables': 'apg_ts_s',
        '_FillValue': np.nan
    }
})


class Calibrations(Coefficients):
    def __init__(self, coeff_file):
        """
        Loads the OPTAA factory calibration coefficients for a unit. Values
        come from either a serialized object created per instrument and
        deployment (calibration coefficients do not change in the middle of a
        deployment), or from the calibration data available from the OOI M2M
        API.
        """
        # assign the inputs
        Coefficients.__init__(self, coeff_file)

    def parse_m2m_cals(self, serial_number, cals, cal_idx):
        """
        Parse the calibration data from the M2M object store. The calibration
        data is stored as an unsorted list of dictionaries. The cal_idx is
        used to identify the calibration data of interest in each dictionary.

        :param serial_number: instrument serial number
        :param cals: list of calibration dictionaries
        :param cal_idx: dictionary of calibration indices
        :return: dictionary of calibration coefficients
        """
        # create the device file dictionary and assign values
        coeffs = {}

        # parse the calibration data
        for cal in cals:
            # beam attenuation and absorption channel clear water offsets
            if cal['name'] == 'CC_acwo':
                coeffs['a_offsets'] = np.array(cal['calData'][cal_idx['CC_acwo']]['value'])
            if cal['name'] == 'CC_ccwo':
                coeffs['c_offsets'] = np.array(cal['calData'][cal_idx['CC_ccwo']]['value'])
            # beam attenuation and absorption channel wavelengths
            if cal['name'] == 'CC_awlngth':
                coeffs['a_wavelengths'] = np.array(cal['calData'][cal_idx['CC_awlngth']]['value'])
            if cal['name'] == 'CC_cwlngth':
                coeffs['c_wavelengths'] = np.array(cal['calData'][cal_idx['CC_cwlngth']]['value'])
            # internal temperature compensation values
            if cal['name'] == 'CC_tbins':
                coeffs['temp_bins'] = np.array(cal['calData'][cal_idx['CC_tbins']]['value'])
            # temperature of calibration water
            if cal['name'] == 'CC_tcal':
                coeffs['temp_calibration'] = cal['calData'][cal_idx['CC_tcal']]['value']
            # temperature compensation values as f(wavelength, temperature) for the attenuation and absorption channels
            if cal['name'] == 'CC_tcarray':
                coeffs['tc_array'] = np.array(cal['calData'][cal_idx['CC_tcarray']]['value'])
            if cal['name'] == 'CC_taarray':
                coeffs['ta_array'] = np.array(cal['calData'][cal_idx['CC_taarray']]['value'])

        # number of wavelengths
        coeffs['num_wavelengths'] = len(coeffs['a_wavelengths'])
        # number of internal temperature compensation bins
        coeffs['num_temp_bins'] = len(coeffs['temp_bins'])
        # pressure coefficients, set to 0 since not included in the CI csv files
        coeffs['pressure_coeff'] = [0, 0]

        # determine the grating index
        awlngths = copy(coeffs['a_wavelengths'])
        awlngths[(awlngths < 545) | (awlngths > 605)] = np.nan
        cwlngths = copy(coeffs['c_wavelengths'])
        cwlngths[(cwlngths < 545) | (cwlngths > 605)] = np.nan
        grate_index = np.nanargmin(np.diff(awlngths) + np.diff(cwlngths))
        coeffs['grate_index'] = grate_index

        # serial number, stripping off all but the numbers
        coeffs['serial_number'] = int(re.sub('[^0-9]', '', serial_number))

        # save the resulting dictionary
        self.coeffs = coeffs


def load_cal_coefficients(cal_file, uid, start_time):
    """
    Load the calibration coefficients for the instrument and deployment from
    the OOI M2M system or from a local file. If the local file does not exist,
    the calibration coefficients will be downloaded from the OOI M2M system and
    saved to the local file for future use.

    :param cal_file: path to the local calibration file
    :param uid: instrument unique identifier (UID)
    :param start_time: deployment start time in seconds since 1970-01-01
    :return: calibration coefficients dictionary
    """
    # load the instrument calibration data
    dev = Calibrations(cal_file)  # initialize calibration class

    # check for the source of calibration coeffs and load accordingly
    if os.path.isfile(cal_file):
        # we always want to use this file if it already exists
        dev.load_coeffs()
    else:
        # load from the OOI M2M system and create a list of calibration events relative to the deployment start date
        cals = get_calibrations_by_uid(uid)
        cal_idx = {}
        for cal in cals['calibration']:
            # for each instance of a calibration event for this instrument...
            tdiff = []
            for data in cal['calData']:
                # calculate the time difference between the start of the deployment and the calibration event
                tdiff.append(data['eventStartTime'] / 1000 - start_time)

            # find the calibration event closest to the start of the deployment
            cal_idx[cal['name']] = np.argmin(np.abs(tdiff))

        # load the calibration coefficients
        dev.parse_m2m_cals(cals['serialNumber'], cals['calibration'], cal_idx)

        # save the calibration coefficients
        dev.save_coeffs()

    return dev


@dask.delayed
def _pd_calc(ref, sig, offset, tintrn, tbins, tarray):
    """
    Convert the raw reference and signal measurements to scientific units. Uses
    a simplified version of the opt_pd_calc function from the pyseas library
    (a fork of the OOI ion_functions code converted to Python 3) to take
    advantage of numpy arrays and the ability to vectorize some of the
    calculations.

    :param ref: raw reference light measurements (OPTCREF_L0 or OPTAREF_L0, as
        appropriate) [counts]
    :param sig: raw signal light measurements (OPTCSIG_L0 or OPTASIG_L0, as
            appropriate) [counts]
    :param offset: 'a' or 'c' (as appropriate) clear water offsets from the
        AC-S device file [m-1]
    :param tintrn: internal instrument temperature [deg_C]
    :param tbins: instrument specific internal temperature calibration bin
        values from AC-S device file [deg_C]
    :param tarray: instrument, wavelength and channel ('c' or 'a') specific
        internal temperature calibration correction coefficients from AC-S
        device file [m-1]
    :return: uncorrected beam attenuation/optical absorption coefficients [m-1]
    """
    # create a linear temperature correction factor based on the internal instrument temperature
    # find the temperature bins corresponding to the values bracketing the internal temperature.
    t0 = tbins[tbins - tintrn < 0][-1]  # set first bracketing temperature
    t1 = tbins[tintrn - tbins < 0][0]   # set second bracketing temperature

    # use the temperature bins to select the calibration coefficients bracketing the internal
    tbins = list(tbins)
    dt0 = tarray[:, tbins.index(t0)]
    dt1 = tarray[:, tbins.index(t1)]

    # Calculate the linear temperature correction.
    tcorr = dt0 + ((tintrn - t0) / (t1 - t0)) * (dt1 - dt0)

    # convert the raw data to the uncorrected spectra
    pd = (offset - (1. / 0.25) * np.log(sig / ref)) - tcorr
    return pd


@dask.delayed(nout=2)
def _holo_grater(wlngths, spectra, index):
    """
    Derived from the Matlab HoloGrater function in Jesse Bausell's
    acsPROCESS_INTERACTIVE toolbox (link below) used in preparing
    AC-S data for NASA's SeaBASS submission process. From the
    original source:

    "This function performs the holographic grating correction for raw
    AC-S spectra. For each individual spectrum it calculates expected
    absorption/attenuation at the lowest wavelength of the second grating
    (upper wavelengths) using matlab's spline function. It then subtracts
    this value from the observed absorption/attenuation creating an offset."

    This function utilizes the SciPy CubicSpline function to accomplish the
    same correction.

    For original code see https://github.com/JesseBausell/acsPROCESS_INTERACTIVE

    :param wlngths: absorption or attenuation channel wavelengths [nm]
    :param spectra: absorption or attenuation spectra [m-1]
    :param index: index of the second holographic grating
    :return: corrected spectra and offset
    """
    # Interpolate between holographic gratings and calculate the offset
    spl = CubicSpline(wlngths[index - 2:index + 1], spectra[index - 2:index + 1])
    interpolation = spl(wlngths[index + 1], extrapolate=True)

    # calculate the offset as the difference between the observed and expected absorption/attenuation
    offset = interpolation - spectra[index + 1]

    # use the offset to correct the second holographic grating
    spectra[index + 1:] = spectra[index + 1:] + offset

    return spectra, offset


def apply_dev(optaa, coeffs):
    """
    Processes the raw data contained in the optaa dictionary and applies the
    factory calibration coefficients contained in the coeffs dictionary to
    convert the data into initial science units. Processing includes correcting
    for the holographic grating offset common to AC-S instruments.

    :param optaa: xarray dataset with the raw absorption and beam attenuation
        measurements.
    :param coeffs: Factory calibration coefficients in a dictionary structure

    :return optaa: xarray dataset with the raw absorption and beam attenuation
        measurements converted into particulate and beam attenuation values
        with the factory pure water calibration values subtracted.
    """
    # calculate the L1 OPTAA data products (uncorrected beam attenuation and absorbance) for particulate
    # and dissolved organic matter with pure water removed.
    nrows = optaa['a_reference'].shape[0]
    pg = []
    a_ref = optaa['a_reference'].values
    a_sig = optaa['a_signal'].values
    c_ref = optaa['c_reference'].values
    c_sig = optaa['c_signal'].values
    temp_internal = optaa['internal_temp'].values
    with tqdm(total=nrows, desc="Collecting data prior to computing the L1 data products") as pbar:
        for i in range(nrows):
            pg += _pd_calc(a_ref[i, :], a_sig[i, :], coeffs['a_offsets'], temp_internal[i],
                           coeffs['temp_bins'], coeffs['ta_array']),
            pg += _pd_calc(c_ref[i, :], c_sig[i, :], coeffs['c_offsets'], temp_internal[i],
                           coeffs['temp_bins'], coeffs['tc_array']),
            pbar.update()

    with ProgressBar():
        print("Calculating the L1 data products for the absorption and attenuation channels")
        pg = [*dask.compute(*pg, scheduler='threads')]

    # create data arrays of the L1 data products
    apg = np.array(pg[0::2])
    m = ~np.isfinite(apg)
    apg[m] = np.nan
    cpg = np.array(pg[1::2])
    m = ~np.isfinite(cpg)
    cpg[m] = np.nan

    # correct the spectral "jump" often observed between the two halves of the linear variable filter
    if coeffs['grate_index']:
        jumps = []
        with tqdm(total=nrows, desc="Collecting data to adjust for spectral jumps between filter halves") as pbar:
            for i in range(nrows):
                spectra, offset = _holo_grater(coeffs['a_wavelengths'], apg[i, :], coeffs['grate_index'])
                jumps += spectra,
                jumps += offset,
                spectra, offset = _holo_grater(coeffs['c_wavelengths'], cpg[i, :], coeffs['grate_index'])
                jumps += spectra,
                jumps += offset,
                pbar.update()

        # put it all back together, adding the jump offsets to the data set
        with ProgressBar():
            print("Adjusting the spectra for the jump often observed between filter halves")
            jumps = [*dask.compute(*jumps, scheduler='threads')]

        optaa['a_jump_offsets'] = ('time', np.array(jumps[1::4]))
        optaa['c_jump_offsets'] = ('time', np.array(jumps[3::4]))

    # return the L1 data with the factory calibrations applied and the spectral jump corrected (if available)
    optaa['apg'] = (('time', 'wavelength_number'), np.array(jumps[0::4]))
    optaa['cpg'] = (('time', 'wavelength_number'), np.array(jumps[2::4]))
    return optaa


def tempsal_corr(channel, pd, wlngth, tcal, temperature, salinity):
    """
    Apply temperature and salinity corrections to the converted absorption
    and attenuation data. Uses a simplified version of the opt_tempsal_corr
    function from the pyseas library (fork of the OOI ion_functions code
    converted to Python 3) to take advantage of numpy arrays and the
    ability to "vectorize" some of the calculations.

    :param channel: string ('a' or 'c') indicating either the absorption or
        attenuation channel is being corrected
    :param pd: array of converted absorption or attenuation data
    :param wlngth: absorption or attenuation channel wavelengths from the
        calibration coefficients
    :param tcal: temperature of the pure water used in the calibrations
    :param temperature: in-situ temperature, ideally from a co-located CTD
    :param salinity: in-situ salinity, ideally from a co-located CTD
    :return: temperature and salinity corrected data
    """
    # create the temperature and salinity correction arrays for each wavelength
    cor_coeffs = np.array([tscor[ii] for ii in wlngth])
    nrows = len(temperature)

    temp_corr = np.tile(cor_coeffs[:, 0], [nrows, 1])
    saln_c_corr = np.tile(cor_coeffs[:, 1], [nrows, 1])
    saln_a_corr = np.tile(cor_coeffs[:, 2], [nrows, 1])

    delta_temp = np.atleast_2d(temperature - tcal).T
    salinity = np.atleast_2d(salinity).T

    if channel == 'a':
        pd_ts = pd - delta_temp * temp_corr - salinity * saln_a_corr
    elif channel == 'c':
        pd_ts = pd - delta_temp * temp_corr - salinity * saln_c_corr
    else:
        raise ValueError('Channel must be either "a" or "c"')

    return pd_ts


def apply_tscorr(optaa, coeffs, temperature, salinity):
    """
    Corrects the absorption and beam attenuation data for the absorption
    of seawater as a function of the seawater temperature and salinity (the
    calibration blanking offsets are determined using pure water.)

    If inputs temperature or salinity are not supplied as calling arguments,
    or all of the temperature or salinity values are NaN, then the following
    default values are used.

        temperature: temperature values recorded by the AC-S's external
            thermistor (note, this would not be valid for an AC-S on a
            profiling platform)
        salinity: 34.0 psu

    Otherwise, each of the arguments for temp and salinity should be either a
    scalar, or a 1D array or a row or column vector with the same number of time
    points as 'a' and 'c'.

    :param optaa: xarray dataset with the raw absorption and attenuation data
        converted to absorption and attenuation coefficients
    :param coeffs: Factory calibration coefficients in a dictionary structure
    :param temperature: In-situ seawater temperature, from a co-located CTD
    :param salinity: In-situ seawater salinity, from a co-located CTD

    :return optaa: xarray dataset with the temperature and salinity corrected
        absorbance and attenuation data arrays added.
    """
    # check the temperature and salinity inputs. If they are not supplied, use the
    # external thermistor temperature and a salinity of 34.0 psu
    if temperature is None:
        temperature = optaa['external_temp']
    if salinity is None:
        salinity = np.ones_like(temperature) * 34.0

    # additionally check if all the temperature and salinity values are NaNs. If they are,
    # then use the external thermistor temperature and a salinity of 34.0 psu (will occur
    # if the CTD is not connected).
    if np.all(np.isnan(temperature)):
        temperature = optaa['external_temp']
    if np.all(np.isnan(salinity)):
        salinity = np.ones_like(temperature) * 34.0

    # test if the temperature and salinity are the same size as the absorption and attenuation
    # data. If they are not, then they should be a scalar value, and we can tile them to the
    # correct size.
    if temperature.size != optaa['time'].size:
        temperature = np.tile(temperature, optaa['time'].size)
    if salinity.size != optaa['time'].size:
        salinity = np.tile(salinity, optaa['time'].size)

    # apply the temperature and salinity corrections
    optaa['apg_ts'] = tempsal_corr('a', optaa['apg'], coeffs['a_wavelengths'], coeffs['temp_calibration'],
                                   temperature, salinity)
    optaa['cpg_ts'] = tempsal_corr('c', optaa['cpg'], coeffs['c_wavelengths'], coeffs['temp_calibration'],
                                   temperature, salinity)

    return optaa


def apply_scatcorr(optaa, coeffs):
    """
    Correct the absorbance data for scattering using Method 1, with the
    wavelength closest to 715 nm used as the reference wavelength for the
    scattering correction. This is the simplest method for correcting for
    scattering, but other methods are available. Users are encouraged to
    explore the other methods and determine which is best for their
    application.

    :param optaa: xarray dataset with the temperature and salinity corrected
        absorbance data array that will be corrected for the effects of
        scattering.
    :param coeffs: Factory calibration coefficients in a dictionary structure

    :return optaa: xarray dataset with the method 1 scatter corrected
        absorbance data array added.
    """
    # find the closest wavelength to 715 nm
    reference_wavelength = 715.0
    idx = np.argmin(np.abs(coeffs['a_wavelengths'] - reference_wavelength))

    # use that wavelength as our scatter correction wavelength
    apg_ts = optaa['apg_ts']
    optaa['apg_ts_s'] = apg_ts - apg_ts[:, idx]

    return optaa


def estimate_chl_poc(optaa, coeffs, chl_line_height=0.020):
    """
    Derive estimates of Chlorophyll-a and particulate organic carbon (POC)
    concentrations from the temperature, salinity and scatter corrected
    absorption and beam attenuation data.

    :param optaa: xarray dataset with the scatter corrected absorbance data.
    :param coeffs: Factory calibration coefficients in a dictionary structure
    :param chl_line_height: Extinction coefficient for estimating the
        chlorophyll concentration. This value may vary regionally and/or
        seasonally. A default value of 0.020 is used if one is not entered,
        but users may to adjust this based on cross-comparisons with other
        measures of chlorophyll
    :return optaa: xarray dataset with the estimates for chlorophyll and POC
        concentrations added.
    """
    # use the standard chlorophyll line height estimation with an extinction coefficient of 0.020.
    m650 = np.argmin(np.abs(coeffs['a_wavelengths'] - 650.0))  # find the closest wavelength to 650 nm
    m676 = np.argmin(np.abs(coeffs['a_wavelengths'] - 676.0))  # find the closest wavelength to 676 nm
    m715 = np.argmin(np.abs(coeffs['a_wavelengths'] - 715.0))  # find the closest wavelength to 715 nm
    apg = optaa['apg_ts_s']
    aphi = apg[:, m676] - 39 / 65 * apg[:, m650] - 26 / 65 * apg[:, m715]
    optaa['estimated_chlorophyll'] = aphi / chl_line_height

    # estimate the POC concentration from the attenuation at 660 nm
    m660 = np.argmin(np.abs(coeffs['c_wavelengths'] - 660.0))  # find the closest wavelength to 660 nm
    cpg = optaa['cpg_ts']
    optaa['estimated_poc'] = cpg[:, m660] * 380

    return optaa


def calculate_ratios(optaa):
    """
    Pigment ratios can be calculated to assess the impacts of bio-fouling,
    sensor calibration drift, potential changes in community composition,
    light history or bloom health and age. Calculated ratios are:

    * CDOM Ratio -- ratio of CDOM absorption in the violet portion of the
        spectrum at 412 nm relative to chlorophyll absorption at 440 nm.
        Ratios greater than 1 indicate a preponderance of CDOM absorption
        relative to chlorophyll.
    * Carotenoid Ratio -- ratio of carotenoid absorption in the blue-green
        portion of the spectrum at 490 nm relative to chlorophyll absorption at
        440 nm. A changing carotenoid to chlorophyll ratio may indicate a shift
        in phytoplankton community composition in addition to changes in light
        history or bloom health and age.
    * Phycobilin Ratio -- ratio of phycobilin absorption in the green portion
        of the spectrum at 530 nm relative to chlorophyll absorption at 440 nm.
        Different phytoplankton, notably cyanobacteria, utilize phycobilins as
        accessory light harvesting pigments. An increasing phycobilin to
        chlorophyll ratio may indicate a shift in phytoplankton community
        composition.
    * Q Band Ratio -- the Soret and the Q bands represent the two main
        absorption bands of chlorophyll. The former covers absorption in the
        blue region of the spectrum, while the latter covers absorption in the
        red region. A decrease in the ratio of the intensity of the Soret band
        at 440 nm to that of the Q band at 676 nm may indicate a change in
        phytoplankton community structure. All phytoplankton contain
        chlorophyll 'a' as the primary light harvesting pigment, but green
        algae and dinoflagellates contain chlorophyll 'b' and 'c', respectively,
        which are spectrally redshifted compared to chlorophyll 'a'.

    :param optaa: xarray dataset with the scatter corrected absorbance data.
    :return optaa: xarray dataset with the estimates for chlorophyll and POC
        concentrations added.
    """
    apg = optaa['optical_absorption']
    m412 = np.nanargmin(np.abs(optaa['wavelength_a'].values[0, :] - 412.0))
    m440 = np.nanargmin(np.abs(optaa['wavelength_a'].values[0, :] - 440.0))
    m490 = np.nanargmin(np.abs(optaa['wavelength_a'].values[0, :] - 490.0))
    m530 = np.nanargmin(np.abs(optaa['wavelength_a'].values[0, :] - 530.0))
    m676 = np.nanargmin(np.abs(optaa['wavelength_a'].values[0, :] - 676.0))

    optaa['ratio_cdom'] = apg[:, m412] / apg[:, m440]
    optaa['ratio_carotenoids'] = apg[:, m490] / apg[:, m440]
    optaa['ratio_phycobilins'] = apg[:, m530] / apg[:, m440]
    optaa['ratio_qband'] = apg[:, m676] / apg[:, m440]

    return optaa


def optaa_datalogger(ds, cal_file):
    """
    Takes OPTAA data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-calculate the intermediate products (e.g. absorption and attenuation) and
    add them to the data set.  Finally, add the estimated chlorophyll and POC
    concentrations to the data set.

    Will test the data set to determine if more than one deployment is present.
    If so, will raise an exception with an error message.  AC-S processing
    requires that the data be processed one deployment at a time in order to
    properly assign calibration coefficients.

    :param ds: initial optaa data set downloaded from OOI via the M2M system
    :param cal_file: file name (can include path) to store the calibration
        coefficients
    :return ds: cleaned up data set
    """
    # check to see if there is more than one deployment in the data set
    if len(np.unique(ds['deployment'].values)) > 1:
        raise ValueError('More than one deployment in the data set.  Please structure processing request to process '
                         'one deployment at a time.')

    # drop some of the variables:
    #   internal_timestamp == time, redundant so can remove
    #   pressure_counts == none of the OOI OPTAAs have a pressure sensor
    ds = ds.drop(['internal_timestamp', 'pressure_counts'])

    # check for data from a co-located CTD, if not present create the variables using NaN's as the fill value
    if 'sea_water_temperature' not in ds.variables:
        ds['sea_water_temperature'] = ('time', ds['deployment'].data * np.nan)
        ds['sea_water_practical_salinity'] = ('time', ds['deployment'].data * np.nan)

    # pull out the number of wavelengths and serial number and then drop the variable (part of the metadata)
    num_wavelengths = ds.num_wavelengths.values[0].astype(int)
    serial_number = int(re.sub('[^0-9]', '', ds.attrs['SerialNumber']))
    ds = ds.drop('num_wavelengths')

    # load the calibration coefficients
    uid = ds.attrs['AssetUniqueID']
    start_time = ds['time'][0].values.astype(float) / 10 ** 9
    cal = load_cal_coefficients(cal_file, uid, start_time)

    # check the calibration coefficients against the deployment data
    if cal.coeffs['serial_number'] != serial_number:
        raise Exception('Serial Number mismatch between ac-s data and the device file.')
    if cal.coeffs['num_wavelengths'] != num_wavelengths:
        raise Exception('Number of wavelengths mismatch between ac-s data and the device file.')

    # remove the units from the variable names
    rename = {
        'a_signal_dark_counts': 'a_signal_dark',
        'a_reference_dark_counts': 'a_reference_dark',
        'a_signal_counts': 'a_signal',
        'a_reference_counts': 'a_reference',
        'c_signal_dark_counts': 'c_signal_dark',
        'c_reference_dark_counts': 'c_reference_dark',
        'c_signal_counts': 'c_signal',
        'c_reference_counts': 'c_reference',
        'wavelength': 'wavelength_number'
    }
    ds = ds.rename(rename)

    # Delete the first 60 seconds of the data record per recommendation from the vendor
    ds.elapsed_run_time.values = ds.elapsed_run_time.where(ds.elapsed_run_time / 1000 > 60)
    ds = ds.dropna(dim='time', subset=['elapsed_run_time'])

    # convert internal and external temperature sensors from raw counts to degrees Celsius
    ds['internal_temp'] = opt_internal_temp(ds['internal_temp_raw'])
    ds['external_temp'] = opt_external_temp(ds['external_temp_raw'])

    # calculate the median of the remaining data per burst measurement
    print('Calculating burst averages...')
    burst = ds.resample(time='900s', base=3150, loffset='450s', skipna=True).reduce(np.median, dim='time',
                                                                                    keep_attrs=True)
    burst = burst.where(~np.isnan(burst.deployment), drop=True)

    # re-process the raw data in order to create the intermediate variables, correcting for the holographic
    # grating, applying the temperature and salinity corrections and applying a baseline scatter correction
    # to the absorption data. All intermediate processing outputs are added to the data set.
    burst = apply_dev(burst, cal.coeffs)
    burst = apply_tscorr(burst, cal.coeffs, burst.sea_water_temperature, burst.sea_water_practical_salinity)
    burst = apply_scatcorr(burst, cal.coeffs)

    # add the jump offsets as NaN's if the grating index correction was not used
    if 'a_jump_offsets' not in ds.variables:
        ds['a_jump_offsets'] = ('time', ds['deployment'].data * np.nan)
        ds['c_jump_offsets'] = ('time', ds['deployment'].data * np.nan)

    # estimate chlorophyll and POC and calculate select absorption ratios
    burst = estimate_chl_poc(burst, cal.coeffs)
    burst = calculate_ratios(burst)

    # create a xarray dataset of the 2D variables, padding the number of wavelengths to a consistent
    # length of 100 using fill values.
    wavelength_number = np.arange(100).astype(int)  # used as a dimensional variable
    pad = 100 - num_wavelengths
    fill_nan = np.tile(np.ones(pad) * np.nan, (len(burst.time), 1))
    fill_int = np.tile(np.ones(pad) * -9999999, (len(burst.time), 1))

    wavelength_a = np.concatenate([burst.wavelength_a.values, fill_nan], axis=1)
    wavelength_c = np.concatenate([burst.wavelength_c.values, fill_nan], axis=1)

    ac = xr.Dataset({
        'wavelength_a': (['time', 'wavelength_number'], wavelength_a),
        'a_signal': (['time', 'wavelength_number'], np.concatenate([burst.a_signal.astype(int), fill_int], axis=1)),
        'a_reference': (['time', 'wavelength_number'], np.concatenate([burst.a_reference.astype(int), fill_int],
                                                                      axis=1)),
        'optical_absorption': (['time', 'wavelength_number'], np.concatenate([burst.optical_absorption, fill_nan],
                                                                             axis=1)),
        'apg': (['time', 'wavelength_number'], np.concatenate([burst.apg, fill_nan], axis=1)),
        'apg_ts': (['time', 'wavelength_number'], np.concatenate([burst.apg_ts, fill_nan], axis=1)),
        'apg_ts_s': (['time', 'wavelength_number'], np.concatenate([burst.apg_ts_s, fill_nan], axis=1)),
        'wavelength_c': (['time', 'wavelength_number'], wavelength_c),
        'c_signal': (['time', 'wavelength_number'], np.concatenate([burst.c_signal.astype(int), fill_int], axis=1)),
        'c_reference': (['time', 'wavelength_number'], np.concatenate([burst.c_reference.astype(int), fill_int],
                                                                      axis=1)),
        'beam_attenuation': (['time', 'wavelength_number'], np.concatenate([burst.beam_attenuation, fill_nan], axis=1)),
        'cpg': (['time', 'wavelength_number'], np.concatenate([burst.cpg, fill_nan], axis=1)),
        'cpg_ts': (['time', 'wavelength_number'], np.concatenate([burst.cpg_ts, fill_nan], axis=1)),
    }, coords={'time': (['time'], burst.time.values), 'wavelength_number': wavelength_number})

    # drop the original 2D variables from the burst data set
    drop = burst.drop(['wavelength_number', 'wavelength_a', 'a_signal', 'a_reference',
                       'optical_absorption', 'apg', 'apg_ts', 'apg_ts_s',
                       'wavelength_c', 'c_signal', 'c_reference',
                       'beam_attenuation', 'cpg', 'cpg_ts'])

    # reset the data type for the 'a' and 'c' signal and reference dark values, and the other raw parameters
    int_arrays = ['a_signal_dark', 'a_reference_dark', 'c_signal_dark', 'c_reference_dark',
                  'internal_temp_raw', 'external_temp_raw', 'deployment']
    for k in drop.variables:
        if k in int_arrays:
            drop[k] = drop[k].astype(int)

    # recombine the two datasets
    optaa = xr.merge([drop, ac])

    # reset the attributes, which the merging drops
    optaa.attrs = burst.attrs
    for v in optaa.variables:
        optaa[v].attrs = burst[v].attrs

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in optaa.variables:
                optaa[key].attrs[atk] = atv

    # add the original variable name as an attribute, if renamed
    for key, value in rename.items():
        optaa[value].attrs['ooinet_variable_name'] = key

    # add the actual number of wavelengths to the dataset as an attribute
    optaa['wavelength_number'].attrs['actual_wavelengths'] = num_wavelengths

    # if the filter index was used to adjust the spectral jumps, add that attribute to the data set
    if cal.coeffs['grate_index']:
        optaa['a_jump_offsets'].attrs['grate_index'] = cal.coeffs['grate_index']
        optaa['c_jump_offsets'].attrs['grate_index'] = cal.coeffs['grate_index']

    return optaa


def optaa_cspp(ds, cal_file):
    """
    Takes OPTAA data recorded by the Coastal Surface-Piercing Profiler (CSPP)
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use.  Additionally,
    re-calculate the intermediate products (e.g. absorption and attenuation) and
    add them to the data set.  Finally, add the estimated chlorophyll and POC
    concentrations to the data set.

    Will test the data set to determine if more than one deployment is present.
    If so, will raise an exception with an error message.  AC-S processing
    requires that the data be processed one deployment at a time in order to
    properly assign calibration coefficients.

    :param ds: initial optaa data set downloaded from OOI via the M2M system
    :param cal_file: file name (can include path) to store the calibration
        coefficients
    :return ds: cleaned up data set
    """
    # check to see if there is more than one deployment in the data set
    if len(np.unique(ds['deployment'].values)) > 1:
        raise ValueError('More than one deployment in the data set.  Please structure processing request to process '
                         'one deployment at a time.')

    # drop some of the variables:
    #   internal_timestamp == time, redundant so can remove
    #   profiler_timestamp == internal_timestamp == time, redundant so can remove
    #   suspect_timestamp = not used
    #   pressure_counts == none of the OOI OPTAAs have a pressure sensor
    ds = ds.drop(['internal_timestamp', 'profiler_timestamp', 'suspect_timestamp', 'pressure_counts'])

    # check for data from a co-located CTD, if not present create the variables using NaN's as the fill value
    if 'sea_water_temperature' not in ds.variables:
        ds['sea_water_temperature'] = ('time', ds['deployment'].data * np.nan)
        ds['sea_water_practical_salinity'] = ('time', ds['deployment'].data * np.nan)

    # pull out the number of wavelengths and serial number and then drop the variable (part of the metadata)
    num_wavelengths = ds.num_wavelengths.values[0].astype(int)
    serial_number = int(re.sub('[^0-9]', '', ds.attrs['SerialNumber']))
    ds = ds.drop('num_wavelengths')

    # load the calibration coefficients
    uid = ds.attrs['AssetUniqueID']
    start_time = ds['time'][0].values.astype(float) / 10 ** 9
    cal = load_cal_coefficients(cal_file, uid, start_time)

    # check the calibration coefficients against the deployment data
    if cal.coeffs['serial_number'] != serial_number:
        raise Exception('Serial Number mismatch between ac-s data and the device file.')
    if cal.coeffs['num_wavelengths'] != num_wavelengths:
        raise Exception('Number of wavelengths mismatch between ac-s data and the device file.')

    # remove the units from the variable name
    rename = {
        'a_signal_dark_counts': 'a_signal_dark',
        'a_reference_dark_counts': 'a_reference_dark',
        'a_signal_counts': 'a_signal',
        'a_reference_counts': 'a_reference',
        'c_signal_dark_counts': 'c_signal_dark',
        'c_reference_dark_counts': 'c_reference_dark',
        'c_signal_counts': 'c_signal',
        'c_reference_counts': 'c_reference',
        'on_seconds': 'elapsed_run_time',
        'int_ctd_pressure': 'sea_water_pressure',
        'wavelength': 'wavelength_number'
    }
    ds = ds.rename(rename)

    # Delete the first 60 seconds of the data record per recommendation from the vendor
    ds.elapsed_run_time.values = ds.elapsed_run_time.where(ds.elapsed_run_time > 60)
    ds = ds.dropna(dim='time', subset=['elapsed_run_time'])

    # convert internal and external temperature sensors from raw counts to degrees Celsius
    ds['internal_temp'] = opt_internal_temp(ds['internal_temp_raw'])
    ds['external_temp'] = opt_external_temp(ds['external_temp_raw'])

    # create a profile variable to uniquely identify profiles within the dataset
    print('Creating and adding a profile variable to the data set.')
    ds = create_profile_id(ds)

    # group the data by profile number and bin the data into 25 cm depth bins (nominal ascent rate of the CSPP)
    vocab = get_vocabulary(ds.attrs['subsite'], ds.attrs['node'], ds.attrs['sensor'])[0]
    site_depth = vocab['maxdepth'] - 2
    profiles = ds.groupby('profile')
    profiles = [profile[1] for profile in profiles]
    partial_binning = partial(bin_profiles, site_depth=site_depth, bin_size=0.25)
    with ProcessPoolExecutor(max_workers=N_CORES) as executor:
        binned = list(tqdm(executor.map(partial_binning, profiles), total=len(profiles),
                           desc='Smoothing and binning each profile into 25 cm depth bins', file=sys.stdout))

    # reset the dataset now using binned profiles
    binned = [i[0] for i in binned if i is not None]
    binned = xr.concat(binned, 'time')
    binned = binned.sortby(['profile', 'time'])

    # confirm dimension order is correct for the wavelength arrays (sometimes the order gets flipped
    # during the binning process)
    binned['wavelength_a'] = binned.wavelength_a.transpose(*['time', 'wavelength_number'])
    binned['wavelength_c'] = binned.wavelength_c.transpose(*['time', 'wavelength_number'])

    # re-process the raw data in order to create the intermediate variables, correcting for the holographic
    # grating, applying the temperature and salinity corrections and applying a baseline scatter correction
    # to the absorption data. All intermediate processing outputs are added to the data set.
    binned = apply_dev(binned, cal.coeffs)
    binned = apply_tscorr(binned, cal.coeffs, binned.sea_water_temperature, binned.sea_water_practical_salinity)
    binned = apply_scatcorr(binned, cal.coeffs)

    # add the jump offsets as NaN's if the grating index correction was not used
    if 'a_jump_offsets' not in ds.variables:
        ds['a_jump_offsets'] = ('time', ds['deployment'].data * np.nan)
        ds['c_jump_offsets'] = ('time', ds['deployment'].data * np.nan)

    # estimate chlorophyll and POC and calculate select absorption ratios
    binned = estimate_chl_poc(binned, cal.coeffs)
    binned = calculate_ratios(binned)

    # create a xarray dataset of the 2D variables, padding the number of wavelengths to a consistent
    # length of 100 using fill values.
    wavelength_number = np.arange(100).astype(int)  # used as a dimensional variable
    pad = 100 - num_wavelengths
    fill_nan = np.tile(np.ones(pad) * np.nan, (len(binned.time), 1))
    fill_int = np.tile(np.ones(pad) * -9999999, (len(binned.time), 1))

    wavelength_a = np.concatenate([binned.wavelength_a.values, fill_nan], axis=1)
    wavelength_c = np.concatenate([binned.wavelength_c.values, fill_nan], axis=1)

    ac = xr.Dataset({
        'wavelength_a': (['time', 'wavelength_number'], wavelength_a),
        'a_signal': (['time', 'wavelength_number'], np.concatenate([binned.a_signal.astype(int), fill_int], axis=1)),
        'a_reference': (['time', 'wavelength_number'], np.concatenate([binned.a_reference.astype(int), fill_int],
                                                                      axis=1)),
        'optical_absorption': (['time', 'wavelength_number'], np.concatenate([binned.optical_absorption, fill_nan],
                                                                             axis=1)),
        'apg': (['time', 'wavelength_number'], np.concatenate([binned.apg, fill_nan], axis=1)),
        'apg_ts': (['time', 'wavelength_number'], np.concatenate([binned.apg_ts, fill_nan], axis=1)),
        'apg_ts_s': (['time', 'wavelength_number'], np.concatenate([binned.apg_ts_s, fill_nan], axis=1)),
        'wavelength_c': (['time', 'wavelength_number'], wavelength_c),
        'c_signal': (['time', 'wavelength_number'], np.concatenate([binned.c_signal.astype(int), fill_int], axis=1)),
        'c_reference': (['time', 'wavelength_number'], np.concatenate([binned.c_reference.astype(int), fill_int],
                                                                      axis=1)),
        'beam_attenuation': (['time', 'wavelength_number'], np.concatenate([binned.beam_attenuation, fill_nan],
                                                                           axis=1)),
        'cpg': (['time', 'wavelength_number'], np.concatenate([binned.cpg, fill_nan], axis=1)),
        'cpg_ts': (['time', 'wavelength_number'], np.concatenate([binned.cpg_ts, fill_nan], axis=1)),
    }, coords={'time': (['time'], binned.time.values), 'wavelength_number': wavelength_number})

    # drop the original 2D variables from the binned data set
    drop = binned.drop(['wavelength_number', 'wavelength_a', 'a_signal', 'a_reference',
                        'optical_absorption', 'apg', 'apg_ts', 'apg_ts_s',
                        'wavelength_c', 'c_signal', 'c_reference',
                        'beam_attenuation', 'cpg', 'cpg_ts'])

    # reset the data type for the 'a' and 'c' signal and reference dark values, and the other raw parameters
    int_arrays = ['a_signal_dark', 'a_reference_dark', 'c_signal_dark', 'c_reference_dark',
                  'internal_temp_raw', 'external_temp_raw', 'deployment', 'profile']
    for k in drop.variables:
        if k in int_arrays:
            drop[k] = drop[k].astype(int)

    # recombine the two datasets
    optaa = xr.merge([drop, ac])

    # reset the attributes, which the merging drops
    optaa.attrs = binned.attrs
    for v in optaa.variables:
        optaa[v].attrs = binned[v].attrs

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in optaa.variables:
                optaa[key].attrs[atk] = atv

    # add the original variable name as an attribute, if renamed
    for key, value in rename.items():
        optaa[value].attrs['ooinet_variable_name'] = key

    # add the actual number of wavelengths to the dataset as an attribute
    optaa['wavelength_number'].attrs['actual_wavelengths'] = num_wavelengths

    # if the filter index was used to adjust the spectral jumps, add that attribute to the data set
    if cal.coeffs['grate_index']:
        optaa['a_jump_offsets'].attrs['grate_index'] = cal.coeffs['grate_index']
        optaa['c_jump_offsets'].attrs['grate_index'] = cal.coeffs['grate_index']

    return optaa


def main(argv=None):
    """
    Command line interface for processing OOI OPTAA NetCDF file(s) from the
    Endurance, Pioneer or Global surface moorings, or the Endurance surface
    piercing profilers. Creates a cleaned and processed xarray dataset of the
    OPTAA data saved to a NetCDF file.
    """
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
        optaa = load_gc_thredds(site, node, sensor, method, stream, ('.*deployment%04d.*OPTAA.*\\.nc$' % deploy))
        cal_file = ('{}.{}.{}.deploy{:02d}.cal_coeffs.json'.format(site, node, sensor, deploy))

        # check to see if we downloaded any data
        if not optaa:
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

        # OPTAA data is different from other instruments. it needs to be processed on a per-deployment basis in order
        # to get the correct number of wavelengths before it can be merged into a single dataset. create a list of
        # all the files that were returned by the M2M request, and determine the deployments that are included in the
        # request
        files = list_files(r['allURLs'][0], '.+OPTAA.+\\.nc$')
        if not files:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

        deployments = np.unique([int(sub.split('/')[3][10:14]) for sub in files])

        # loop through the deployments and download the data for each one
        optaa = []
        cal_file = []
        for deploy in deployments:
            # Valid M2M request, download the data on a per-deployment basis
            data = m2m_collect(r, ('.*deployment%04d.*OPTAA.*\\.nc$' % deploy))
            if data:
                optaa.append(data)
                cal_file.append('{}.{}.{}.deploy{:02d}.cal_coeffs.json'.format(site, node, sensor, deploy))

        # check to see if we downloaded any data (remove empty/none entries from the list)
        if not optaa:
            exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                    stream, start, stop))
            raise SystemExit(exit_text)

    # set up the calibration file path and name(s)
    out_file = os.path.abspath(args.outfile)
    cal_path = os.path.dirname(out_file)
    if not os.path.exists(cal_path):
        os.makedirs(cal_path)

    # clean-up and reorganize the data
    multi = isinstance(optaa, list)
    if node == 'SP001':
        # this OPTAA is part of a CSPP
        if multi:
            for i, ds in enumerate(optaa):
                cfile = os.path.join(cal_path, cal_file[i])
                optaa[i] = optaa_cspp(ds, cfile)
            optaa = xr.concat(optaa, dim='time')
        else:
            cal_file = os.path.join(cal_path, cal_file)
            optaa = optaa_cspp(optaa, cal_file)
    else:
        # this OPTAA is stand-alone on one of the moorings
        if multi:
            for i, ds in enumerate(optaa):
                cfile = os.path.join(cal_path, cal_file[i])
                optaa[i] = optaa_datalogger(ds, cfile)
            optaa = xr.concat(optaa, dim='time')
        else:
            optaa = optaa_datalogger(optaa, cal_file)

    # get the vocabulary information for the site, node, and sensor and update the dataset attributes
    vocab = get_vocabulary(site, node, sensor)[0]
    optaa = optaa.sortby(['deployment', 'time'])
    optaa = update_dataset(optaa, vocab['maxdepth'])

    # save the data to disk
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    optaa.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
