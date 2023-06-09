#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import re

from copy import copy
from functools import partial
from scipy.interpolate import CubicSpline
from tqdm import tqdm

from ooi_data_explorations.calibrations import Coefficients
from ooi_data_explorations.common import get_calibrations_by_uid
from pyseas.data.opt_functions_tscor import tscor


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

    def _compare_names(self, cal_name, data_source):
        if not cal_name:
            cal_name = data_source
        else:
            if cal_name != data_source:
                raise ValueError('Calibration data file name inconsistent, unable to properly parse the '
                                 'calibration data.')

        return cal_name

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
        cal_name = None
        # parse the calibration data
        for cal in cals:
            # beam attenuation and absorption channel clear water offsets
            if cal['name'] == 'CC_acwo':
                coeffs['a_offsets'] = np.array(cal['calData'][cal_idx['CC_acwo']]['value'])
                cal_name = self._compare_names(cal_name, cal['calData'][cal_idx['CC_acwo']]['dataSource'])
            if cal['name'] == 'CC_ccwo':
                coeffs['c_offsets'] = np.array(cal['calData'][cal_idx['CC_ccwo']]['value'])
                cal_name = self._compare_names(cal_name, cal['calData'][cal_idx['CC_ccwo']]['dataSource'])
            # beam attenuation and absorption channel wavelengths
            if cal['name'] == 'CC_awlngth':
                coeffs['a_wavelengths'] = np.array(cal['calData'][cal_idx['CC_awlngth']]['value'])
                cal_name = self._compare_names(cal_name, cal['calData'][cal_idx['CC_awlngth']]['dataSource'])
            if cal['name'] == 'CC_cwlngth':
                coeffs['c_wavelengths'] = np.array(cal['calData'][cal_idx['CC_cwlngth']]['value'])
                cal_name = self._compare_names(cal_name, cal['calData'][cal_idx['CC_cwlngth']]['dataSource'])
            # internal temperature compensation values
            if cal['name'] == 'CC_tbins':
                coeffs['temp_bins'] = np.array(cal['calData'][cal_idx['CC_tbins']]['value'])
                cal_name = self._compare_names(cal_name, cal['calData'][cal_idx['CC_tbins']]['dataSource'])
            # temperature of calibration water
            if cal['name'] == 'CC_tcal':
                coeffs['temp_calibration'] = cal['calData'][cal_idx['CC_tcal']]['value']
                cal_name = self._compare_names(cal_name, cal['calData'][cal_idx['CC_tcal']]['dataSource'])
            # temperature compensation values as f(wavelength, temperature) for the attenuation and absorption channels
            if cal['name'] == 'CC_tcarray':
                coeffs['tc_array'] = np.array(cal['calData'][cal_idx['CC_tcarray']]['value'])
                cal_name = self._compare_names(cal_name, cal['calData'][cal_idx['CC_tcarray']]['dataSource'])
            if cal['name'] == 'CC_taarray':
                coeffs['ta_array'] = np.array(cal['calData'][cal_idx['CC_taarray']]['value'])
                cal_name = self._compare_names(cal_name, cal['calData'][cal_idx['CC_taarray']]['dataSource'])

        # calibration data file name
        coeffs['source_file'] = cal_name
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
                td = (data['eventStartTime'] / 1000) - start_time
                if td <= 0:
                    # valid cals must come before the deployment start date/time
                    tdiff.append(td)
                else:
                    # use a ridiculously large number to avoid selecting this cal event
                    tdiff.append(10**20)

            # find the calibration event closest to the start of the deployment
            cal_idx[cal['name']] = np.argmin(np.abs(tdiff))

        # load the calibration coefficients
        dev.parse_m2m_cals(cals['serialNumber'], cals['calibration'], cal_idx)

        # save the calibration coefficients
        dev.save_coeffs()

    return dev


def convert_raw(ref, sig, tintrn, offset, tarray, tbins):
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
    :param tintrn: internal instrument temperature [deg_C]
    :param offset: 'a' or 'c' (as appropriate) clear water offsets from the
        AC-S device file [m-1]
    :param tarray: instrument, wavelength and channel ('c' or 'a') specific
        internal temperature calibration correction coefficients from AC-S
        device file [m-1]
    :param tbins: instrument specific internal temperature calibration bin
        values from AC-S device file [deg_C]
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
    pg = (offset - (1. / 0.25) * np.log(sig / ref)) - tcorr
    return pg


def holo_grater(wlngths, spectra, index):
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


def pg_calc(reference, signal, internal_temperature, offset, tarray, tbins, grating, wavelengths):
    """
    Combines the convert_raw and holo_grater functions to calculate the L1 data
    products for the AC-S (the uncorrected optical absorption and the beam
    attenuation). This function is intended to be wrapped by the apply_dev
    function to facilitate multiprocessing of the potentially large data arrays
    common to the AC-S.

    :param reference: raw reference light measurements
    :param signal: raw signal light measurements
    :param internal_temperature: internal instrument temperature converted from
        raw counts to degrees Celsius
    :param offset: 'a' or 'c' (as appropriate) clear water offsets from the
        AC-S device file
    :param tarray: instrument, wavelength and channel ('c' or 'a') specific
        internal temperature calibration correction coefficients from AC-S
        device file
    :param tbins: instrument specific internal temperature calibration bin
        values from AC-S device file
    :param grating: index of the for the start of the second half of the
        filter sets
    :param wavelengths: instrument, wavelength and channel ('c' or 'a') specific
        wavelengths from the AC-S device file
    :return: converted optical absorption or beam attenuation measurements
        along with the offset correction for the holographic grating.
    """
    # convert the raw measurements to uncorrected absorption/attenuation values
    pg = convert_raw(reference, signal, internal_temperature, offset, tarray, tbins)
    m = ~np.isfinite(pg)
    pg[m] = np.nan

    # if the grating index is set, correct for the often observed jump at the mid-point of the spectra
    if grating:
        pg, jump = holo_grater(wavelengths, pg, grating)
    else:
        jump = np.nan

    return pg, jump


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
    # pull out the relevant raw data parameters
    nrows = optaa['a_reference'].shape[0]
    a_ref = optaa['a_reference'].values
    a_sig = optaa['a_signal'].values
    c_ref = optaa['c_reference'].values
    c_sig = optaa['c_signal'].values
    temp_internal = optaa['internal_temp'].values

    # create a set of partial functions for concurrent.futures (maps static elements to iterables)
    apg_calc = partial(pg_calc, offset=coeffs['a_offsets'], tarray=coeffs['ta_array'],
                       tbins=coeffs['temp_bins'], grating=coeffs['grate_index'],
                       wavelengths=coeffs['a_wavelengths'])
    cpg_calc = partial(pg_calc, offset=coeffs['c_offsets'], tarray=coeffs['tc_array'],
                       tbins=coeffs['temp_bins'], grating=coeffs['grate_index'],
                       wavelengths=coeffs['c_wavelengths'])

    # apply the partial functions to the data arrays, calculating the L1 data products
    apg = [apg_calc(a_ref[i, :], a_sig[i, :], temp_internal[i]) for i in tqdm(range(nrows),
                                                                              desc='Converting absorption data...')]
    cpg = [cpg_calc(c_ref[i, :], c_sig[i, :], temp_internal[i]) for i in tqdm(range(nrows),
                                                                              desc='Converting attenuation data...')]

    # create data arrays of the L1 data products
    apg, a_jumps = zip(*[row for row in apg])
    apg = np.array(apg)
    m = ~np.isfinite(apg)
    apg[m] = np.nan
    a_jumps = np.array(a_jumps)

    cpg, c_jumps = zip(*[row for row in cpg])
    cpg = np.array(cpg)
    m = ~np.isfinite(cpg)
    cpg[m] = np.nan
    c_jumps = np.array(c_jumps)

    # return the L1 data with the factory calibrations applied and the spectral jump corrected (if available)
    optaa['apg'] = (('time', 'wavelength_number'), apg)
    optaa['cpg'] = (('time', 'wavelength_number'), cpg)
    optaa['a_jump_offsets'] = ('time', a_jumps)
    optaa['c_jump_offsets'] = ('time', c_jumps)
    return optaa


def tempsal_corr(channel, pg, wlngth, tcal, temperature, salinity):
    """
    Apply temperature and salinity corrections to the converted absorption
    and attenuation data. Uses a simplified version of the opt_tempsal_corr
    function from the pyseas library (fork of the OOI ion_functions code
    converted to Python 3) to take advantage of numpy arrays and the
    ability to "vectorize" some of the calculations.

    :param channel: string ('a' or 'c') indicating either the absorption or
        attenuation channel is being corrected
    :param pg: array of converted absorption or attenuation data
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
        pg_ts = pg - delta_temp * temp_corr - salinity * saln_a_corr
    elif channel == 'c':
        pg_ts = pg - delta_temp * temp_corr - salinity * saln_c_corr
    else:
        raise ValueError('Channel must be either "a" or "c"')

    return pg_ts


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
