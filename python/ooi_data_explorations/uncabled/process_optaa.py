#!/usr/bin/env python
# -*- coding: utf-8 -*-
import dateutil.parser as parser
import numpy as np
import os
import pytz
import xarray as xr

from datetime import timedelta

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, load_gc_thredds, \
    list_deployments, get_deployment_dates, get_vocabulary, update_dataset, ENCODINGS

# from cgsn_processing.process.finding_calibrations import find_calibration
# from cgsn_processing.process.proc_optaa import Calibrations, apply_dev, apply_tscorr, apply_scatcorr, \
#     calculate_ratios, estimate_chl_poc

from pyseas.data.opt_functions import opt_internal_temp, opt_external_temp
#from pyseas.data.opt_functions import opt_pd_calc, opt_tempsal_corr

# reset the variable level attributes
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
                    'measures the in-situ seawater termperature.')
    },
    'internal_temp_raw': {
        'long_name': 'Raw Internal Instrument Temperature',
        'units': 'count',
        'comment': ('Raw measurements, reported in counts, from the AC-S internal temperature sensor. This sensor '
                    'measures the internal instrument termperature and is used in converting the raw optical '
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
    'seawater_temperature': {
        'long_name': 'Seawater Temperature',
        'standard_name': 'sea_water_temperature',
        'units': 'degree_Celsius',
        'comment': ('Sea water temperature is the in situ temperature of the sea water. Measurements are from a '
                    'co-located CTD'),
        'data_product_identifier': 'TEMPWAT_L1',
        '_FillValue': np.nan
    },
    'practical_salinity': {
        'long_name': 'Practical Salinity',
        'standard_name': 'sea_water_practical_salinity',
        'units': '1',
        'comment': ('Salinity is generally defined as the concentration of dissolved salt in a parcel of sea water. '
                    'Practical Salinity is a more specific unitless quantity calculated from the conductivity of '
                    'sea water and adjusted for temperature and pressure. It is approximately equivalent to Absolute '
                    'Salinity (the mass fraction of dissolved salt in sea water), but they are not interchangeable. '
                    'Measurements are from a co-located CTD.'),
        'data_product_identifier': 'PRACSAL_L2',
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
    'apd': {
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
    'apd_ts': {
        'long_name': 'Particulate and Dissolved Absorbance with TS Correction',
        'units': 'm-1',
        'comment': ('The optical absorption coefficient corrected for the effects of temperature and salinity. '
                    'Utilizes data from a co-located CTD for the temperaure and salinity, if available. If no '
                    'co-located CTD data is available, will assume a constant salinity of 33 psu and will use '
                    'the OPTAA''s external temperature sensor.'),
        'ancillary_variables': 'wavelength_a ctd_temperature ctd_salinity external_temp apd',
        '_FillValue': np.nan
    },
    'apd_ts_s': {
        'long_name': 'Particulate and Dissolved Absorbance with TS and Scatter Correction',
        'standard_name': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water',
        'units': 'm-1',
        'comment': ('The optical absorption coefficient corrected for scattering after correcting for temperature and '
                    'salinity. Utilizes Method 1 (baseline correction) by subtracting the absorption at 715 nm from '
                    'all values; assumes scattering is flat across wavelengths.'),
        'data_product_identifier': 'OPTABSN_L2',
        'ancillary_variables': 'wavelength_a apd_ts',
        '_FillValue': np.nan
    },
    'cpd': {
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
    'cpd_ts': {
        'long_name': 'Particulate and Dissolved Attenuation with TS Correction',
        'units': 'm-1',
        'comment': ('The optical beam attenuation coefficient corrected for the effects of temperature and salinity. '
                    'Utilizes data from a co-located CTD for the temperaure and salinity, if available. If no '
                    'co-located CTD data is available, will assume a constant salinity of 33 psu and will use '
                    'the OPTAA''s external temperature sensor.'),
        'data_product_identifier': 'OPTATTN_L2',
        'ancillary_variables': 'wavelength_c ctd_temperature ctd_salinity external_temp cpd',
        '_FillValue': np.nan
    },
    'estimated_chlorophyll': {
        'long_name': 'Estimated Chlorophyll Concentration',
        'standard_name': 'mass_concentration_of_chlorophyll_in_sea_water',
        'units': 'ug L-1',
        'comment': ('Uses the absorption line height at 676 nm, above a linear background between 650 and 715 nm with '
                    'a chlorophyll specific absorption of 0.020 L/ug/m, to estimate the concentration of chlorophyll. '
                    'This method has been shown to be significantly related to extracted chlorophyll concentrations '
                    'and is robust in response to mild to moderate biofouling.'),
        'ancillary_variables': 'apd_ts_s',
        '_FillValue': np.nan
    },
    'estimated_poc': {
        'long_name': 'Estimated POC Concentration',
        'standard_name': 'mass_concentration_of_organic_detritus_expressed_as_carbon_in_sea_water',
        'units': 'ug L-1',
        'comment': ('Uses the particulate beam attenuation coefficient at 660 nm and a coefficient of 380 ug/L/m. This '
                    'calculation is not robust in response to biofouling and is expected to breakdown as biofouling '
                    'begins to dominate the signal.'),
        'ancillary_variables': 'cpd_ts',
        '_FillValue': np.nan
    },
    'ratio_cdom': {
        'long_name': 'CDOM to Chlorophyll Absorbance Ratio',
        'units': '1',
        'comment': ('Ratio of CDOM absorption in the violet portion of the spectrum at 412 nm relative to '
                    'chlorophyll absorption at 440 nm. Ratios greater than 1 indicate a preponderance of CDOM '
                    'absorption relative to chlorophyll.'),
        'ancillary_variables': 'apd_ts_s',
        '_FillValue': np.nan
    },
    'ratio_carotenoids': {
        'long_name': 'Carotenoid to Chlorophyll Absorbance Ratio',
        'units': '1',
        'comment': ('Ratio of carotenoid absorption in the blue-green portion of the spectrum at 490 nm relative to '
                    'chlorophyll absorption at 440 nm. A changing carotenoid to chlorophyll ratio may indicate a shift '
                    'in phytoplankton community composition in addition to changes in light history or bloom health '
                    'and age.'),
        'ancillary_variables': 'apd_ts_s',
        '_FillValue': np.nan
    },
    'ratio_phycobilins': {
        'long_name': 'Phycobilins to Chlorophyll Absorbance Ratio',
        'units': '1',
        'comment': ('Ratio of phycobilin absorption in the green portion of the spectrum at 530 nm relative to '
                    'chlorophyll absorption at 440 nm. Different phytoplankton, notably cyanobacteria, utilize '
                    'phycobilins as accessory light harvesting pigments. An increasing phycobilin to chlorophyll ratio '
                    'may indicate a shift in phytoplankton community composition.'),
        'ancillary_variables': 'apd_ts_s',
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
        'ancillary_variables': 'apd_ts_s',
        '_FillValue': np.nan
    }
})


def estimate_chl_poc(optaa):
    """
    The OPTAA data can be used to create estimates of the chlorophyll and
    particulate organic carbon (POC) concentration. These can be compared to
    data from a co-located fluorometer to help validate the performance of
    the sensors.

    Derive estimates of Chlorophyll-a and particulate organic carbon (POC)
    concentrations from the temperature, salinity and scatter corrected
    absorption and beam attenuation data.

    :param optaa: xarray dataset with the scatter corrected absorbance data.
    :return optaa: xarray dataset with the estimates for chlorophyll and POC
        concentrations added.
    """
    # use the standard chlorophyll line height estimation with an extinction coefficient of 0.020.
    m676 = np.nanargmin(np.abs(optaa['wavelength_a'].values[0, :] - 676.0))
    m650 = np.nanargmin(np.abs(optaa['wavelength_a'].values[0, :] - 650.0))
    m715 = np.nanargmin(np.abs(optaa['wavelength_a'].values[0, :] - 715.0))
    apg = optaa['optical_absorption']
    aphi = apg[:, m676] - 39/65 * apg[:, m650] - 26/65 * apg[:, m715]
    optaa['estimated_chlorophyll'] = aphi / 0.020

    # estimate the POC concentration from the attenuation at 660 nm
    m660 = np.nanargmin(np.abs(optaa['wavelength_c'].values[0, :] - 660.0))
    cpg = optaa['beam_attenuation']
    optaa['estimated_poc'] = cpg[:, m660] * 380

    return optaa


def calculate_ratios(optaa):
    """
    Pigment ratios can be calculated to assess the impacts of biofouling,
    sensor calibration drift, potential changes in community composition,
    light history or bloom health and age. Calculated ratios are:

    * CDOM Ratio -- ratio of CDOM absorption in the violet portion of the
        spectrum at 412 nm relative to chlorophyll
    * absorption at 440 nm. Ratios greater than 1 indicate a preponderance of
        CDOM absorption relative to chlorophyll.
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
        chlorophyll a as the primary light harvesting pigment, but green algae
        and dinoflagellates contain chlorophyll b and c, respectively, which
        are spectrally redshifted compared to chlorophyll a.

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


def adjusted_dates(site, node, sensor, deploy):
    """
    Due to a bug in the way the system is currently assigning the calibration
    coefficients, we need to bound the data requests for each deployment
    based on the start and end times of the neighboring deployments. We cannot
    take advantage of the overlapping deployments as that may result the wrong
    calibration coefficients being applied to the record.

    :param site: Site designator, extracted from the first part of the
        reference designator
    :param node: Node designator, extracted from the second part of the
        reference designator
    :param sensor: Sensor designator, extracted from the third and fourth part
        of the reference designator
    :param deploy: sensor deployment number
    :return start: Deployment start date, adjusted to ensure no overlap with
        the previous, if any, deployment
    :return stop: Deployment stop date, adjusted to ensure no overlap with
        the following, if any, deployment
    """
    deployments = list_deployments(site, node, sensor)
    if deploy == deployments[0]:
        # First deployment, use the deployment start date and adjust the stop date
        start, _ = get_deployment_dates(site, node, sensor, deploy)

        stop = parser.parse(get_deployment_dates(site, node, sensor, deploy + 1)[0])
        stop = stop.astimezone(pytz.utc) - timedelta(hours=1)
        stop = stop.strftime('%Y-%m-%dT%H:%M:%S.000Z')
    elif deploy == deployments[-1]:
        # Last deployment, use the deployment end date and adjust the start date
        start = parser.parse(get_deployment_dates(site, node, sensor, deploy - 1)[1])
        start = start.astimezone(pytz.utc) + timedelta(hours=1)
        start = start.strftime('%Y-%m-%dT%H:%M:%S.000Z')

        _, stop = get_deployment_dates(site, node, sensor, deploy)
    else:
        # middle deployments, adjust the start and stop dates
        start = parser.parse(get_deployment_dates(site, node, sensor, deploy - 1)[1])
        start = start.astimezone(pytz.utc) + timedelta(hours=1)
        start = start.strftime('%Y-%m-%dT%H:%M:%S.000Z')

        stop = parser.parse(get_deployment_dates(site, node, sensor, deploy + 1)[0])
        stop = stop.astimezone(pytz.utc) - timedelta(hours=1)
        stop = stop.strftime('%Y-%m-%dT%H:%M:%S.000Z')

    return start, stop


def optaa_datalogger(ds):
    """
    Takes optaa data recorded by the data loggers used in the CGSN/EA moorings
    and cleans up the data set to make it more user-friendly.  Primary task is
    renaming parameters and dropping some that are of limited use. Additionally,
    re-organize some of the variables to permit better assessments of the data.

    :param ds: initial optaa data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    # drop some of the variables:
    #   internal_timestamp == time, redundant so can remove
    #   pressure_counts == none of the OOI OPTAAs have a pressure sensor
    ds = ds.drop(['internal_timestamp', 'pressure_counts'])

    # check for data from a co-located CTD, if not present create the variables using NaN as the fill value
    if 'temp' not in ds.variables:
        ds['temp'] = ('time', ds['deployment'] * np.nan)
        ds['practical_salinity'] = ('time', ds['deployment'] * np.nan)

    # pull out the number of wavelengths and then drop the variable (will add to the metadata)
    num_wavelengths = ds.num_wavelengths.values[0].astype(int)
    ds = ds.drop('num_wavelengths')

    # remove the units from the variable name and rename temp to seawater_temperature
    rename = {
        'a_signal_dark_counts': 'a_signal_dark',
        'a_reference_dark_counts': 'a_reference_dark',
        'a_signal_counts': 'a_signal',
        'a_reference_counts': 'a_reference',
        'c_signal_dark_counts': 'c_signal_dark',
        'c_reference_dark_counts': 'c_reference_dark',
        'c_signal_counts': 'c_signal',
        'c_reference_counts': 'c_reference',
        'temp': 'seawater_temperature',
        'wavelength': 'wavelength_number'
    }
    ds = ds.rename(rename)

    # Delete the first 60 seconds of the data record per recommendation from the vendor
    ds.elapsed_run_time.values = ds.elapsed_run_time.where(ds.elapsed_run_time / 1000 > 60)
    ds = ds.dropna(dim='time', subset=['elapsed_run_time'])

    # calculate the median of the remaining data per burst measurement
    burst = ds.resample(time='900s', base=3150, loffset='450s', skipna=True).reduce(np.median, dim='time',
                                                                                    keep_attrs=True)
    burst = burst.where(~np.isnan(burst.deployment), drop=True)
    burst = burst.compute()

    # create an xarray dataset of the 2D variables, padding the number of wavelengths to a consistent
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
        'wavelength_c': (['time', 'wavelength_number'], wavelength_c),
        'c_signal': (['time', 'wavelength_number'], np.concatenate([burst.c_signal.astype(int), fill_int], axis=1)),
        'c_reference': (['time', 'wavelength_number'], np.concatenate([burst.c_reference.astype(int), fill_int],
                                                                      axis=1)),
        'beam_attenuation': (['time', 'wavelength_number'], np.concatenate([burst.beam_attenuation, fill_nan], axis=1))
    }, coords={'time': (['time'], burst.time.values), 'wavelength_number': wavelength_number})

    # drop the original 2D variables from the burst data set
    drop = burst.drop(['wavelength_number', 'wavelength_a', 'a_signal', 'a_reference', 'optical_absorption',
                        'wavelength_c', 'c_signal', 'c_reference', 'beam_attenuation'])

    # reset the data type for the 'a' and 'c' signal and reference dark values, and the other raw parameters
    int_arrays = ['a_signal_dark', 'a_reference_dark', 'c_signal_dark', 'c_reference_dark',
                  'internal_temp_raw', 'external_temp_raw', 'num_wavelengths', 'deployment']
    for k in drop.variables:
        if k in int_arrays:
            drop[k] = drop[k].astype(int)

    # recombine the two datasets
    optaa = xr.merge([drop, ac])

    # reset the attributes, which the merging drops
    optaa.attrs = ds.attrs
    for v in optaa.variables:
        optaa[v].attrs = ds[v].attrs

    # convert internal and external temperature sensors from raw counts to degrees Celsius
    optaa['internal_temp'] = opt_internal_temp(optaa['internal_temp_raw'])
    optaa['external_temp'] = opt_external_temp(optaa['external_temp_raw'])

    # TODO: recalculate intermediate absorption and attenuation products and add to the dataset

    # Calculate the chlorophyll and particulate organic carbon (POC) concentrations and key pigment ratios
    optaa = estimate_chl_poc(optaa)
    optaa = calculate_ratios(optaa)

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

    return optaa


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

    # if we are specifying a deployment number, then create the start and stop dates
    if deploy:
        start, stop = adjusted_dates(site, node, sensor, deploy)

    # Request the data for download from OOINet via the M2M API using the specified dates
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    if not r:
        exit_text = ('Request failed for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                              stream, start, stop))
        raise SystemExit(exit_text)

    # Valid M2M request, start downloading the data
    optaa = m2m_collect(r, '.*OPTAA.*\\.nc$')

    # check to see if we downloaded any data
    if not optaa:
        exit_text = ('Data unavailable for %s-%s-%s, %s, %s, from %s to %s.' % (site, node, sensor, method,
                                                                                stream, start, stop))
        raise SystemExit(exit_text)

    optaa = optaa_datalogger(optaa)
    vocab = get_vocabulary(site, node, sensor)[0]
    optaa = update_dataset(optaa, vocab['maxdepth'])

    # save the data to disk
    out_file = os.path.abspath(args.outfile)
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    optaa.to_netcdf(out_file, mode='w', format='NETCDF4', engine='h5netcdf', encoding=ENCODINGS)


if __name__ == '__main__':
    main()
