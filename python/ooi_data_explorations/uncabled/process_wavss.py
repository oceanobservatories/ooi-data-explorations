#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os

from ooi_data_explorations.common import inputs, m2m_collect, m2m_request, load_gc_thredds, get_deployment_dates, \
    get_vocabulary, update_dataset, dt64_epoch, ENCODINGS
from ooi_data_explorations.qartod.qc_processing import parse_qc

ATTRS = {
    'average_wave_height': {
        'long_name': 'Mean wave height',
        'standard_name': 'sea_surface_wave_mean_height',
        'units': 'm',
        'comment': ('Wave height is defined as the vertical distance from a wave trough to the following wave crest. '
                    'The mean wave height is the mean trough to crest distance measured during the observation period. '
                    'This is calculated from the average zero down-crossing wave height'),
        'data_product_identifier': 'WAVSTAT-HAVG_L0',
    },
    'max_wave_height': {
        'long_name': 'Maximum wave height',
        'standard_name': 'sea_surface_wave_maximum_height',
        'units': 'm',
        'comment': ('Wave height is defined as the vertical distance from a wave trough to the following wave crest. '
                    'The maximum wave height is the greatest trough to crest distance measured during the observation period. '
                    'This is calculated from the maximum zero down-crossing wave height.'),
        'data_product_identifier': 'WAVSTAT-HMAX_L0',
    },
    'mean_direction': {
        'long_name': 'Mean Direction of Wave Field',
        'standard_name': 'sea_surface_wave_to_direction',
        'units': 'degrees',
        'comment': ('Mean direction of wave field against true north. The phrase "to_direction" is used in the '
                    'construction X_to_direction and indicates the direction towards which the velocity vector of '
                    ' X is headed. The direction is a bearing in the usual geographical sense, measured positive '
                    ' clockwise from due north.'),
        'data_product_identifier': 'WAVSTAT-D_L0',
    },
    'mean_spectral_period': {
        'long_name': 'Mean Wave Period Computed from Second Frequency Moment',
        'standard_name': 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        'units': 's',
        'comment': ('The wave directional spectrum can be written as a five dimensional function S(t,x,y,f,theta) '
                    'where t is time, x and y are horizontal coordinates (such as longitude and latitude), f is '
                    'frequency and theta is direction. S can be integrated over direction to give '
                    'S1 = integral(S dtheta). Frequency moments, M(n) of S1 can then be calculated as follows: '
                    'M(n) = integral(S1 f^n df), where f^n is f to the power of n. '
                    'The second wave period, T(m2) is calculated as the square root of the ratio M(0)/M(2).'),
    },
    'mean_spread': {
        'long_name': 'Mean Directional Spread of the Wave Field',
        'standard_name': 'sea_surface_wave_directional_spread',
        'units': 'degrees',
        'comment': ('Overall directional spreading width in degrees obtained by averaging the spreading '
                    'width sigma theta, σθ, over all frequencies with weighting function S(f). σθ is calculated '
                    'by the KVH method.'),
    },
    'mean_wave_period': {
        'long_name': 'Mean Wave Period',
        'standard_name': 'sea_surface_wave_mean_period',
        'units': 's',
        'comment': ('Wave period is the interval of time between repeated features on the waveform such as crests, '
                    'troughs or upward passes through the mean level. Wave mean period is the mean period measured '
                    'over the observation duration. Calculated as the average zero down-crossing wave period. '),
        'data_product_identifier': 'WAVSTAT-TAVG_L0'
    },
    'peak_wave_period': {
        'long_name': 'Peak Wave Period',
        'standard_name': 'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'units': 's',
        'comment': ('Wave period is the interval of time between repeated features on the waveform such as crests, '
                    'troughs or upward passes through the mean level. The sea_surface_wave_period_at_variance_spectral_density_maximum, '
                    'sometimes called peak wave period, is the period of the most energetic waves in the total wave '
                    'spectrum at a specific location.'),
        'data_product_identifier': 'WAVSTAT-TP_L0',
    },
    'significant_period': {
        'long_name': 'Significant Wave Period',
        'standard_name': 'sea_surface_wave_significant_period',
        'units': 's',
        'comment': ('Wave period is the interval of time between repeated features on the waveform such as crests, '
                    'troughs or upward passes through the mean level. The significant wave period is the mean period '
                    'of the highest one-third of waves. Calculated from the significant zero down-crossing waves.'),
        'data_product_identifier': 'WAVSTAT-TSIG_L0',
    },
    'significant_wave_height': {
        'long_name': 'Significant Wave Height',
        'standard_name': 'sea_surface_wave_significant_height',
        'units': 'm',
        'comment': ('Wave height is defined as the vertical distance from a wave trough to the following wave crest. '
                    'The significant wave height is the mean trough to crest distance measured during the observation '
                    'period of the highest one-third of waves. Calculated from the zero down-crossing significant wave '
                    'height.'),
        'data_product_identifier': 'WAVSTAT-HSIG_L0',
    },
    'wave_height_10': {
        'long_name': 'Mean Wave Height of Highest Tenth',
        'standard_name': 'sea_surface_wave_mean_height_of_highest_tenth',
        'units': 'm',
        'comment': ('Wave height is defined as the vertical distance from a wave trough to the following '
                    'wave crest. The height of the highest tenth is defined as the mean of the highest '
                    'ten per cent of trough to crest distances measured during the observation period.'),
        'data_product_identifier': 'WAVSTAT-H10_L0',
    },
    'wave_height_hmo': {
        'long_name': 'Significant Wave Height from Spectral Moment 0',
        'standard_name': 'sea_surface_wave_significant_height_from_variance_spectral_density',
        'units': 'm',
        'comment': ('Wave height is defined as the vertical distance from a wave trough to the following '
                    'wave crest. The significant wave height (hm0) is the mean wave height of the highest '
                    'one-third of waves as estimated from the zeroth-spectral moment m0, where '
                    'hm0 = 4*sqrt(m0), and m0 is the intregral of the S(f)*df with f = F1 to F2 in Hz'),
        'data_product_identifier': 'WAVSTAT-HMO_L0',
    },
    'wave_period_10': {
        'long_name': 'Mean Wave Period of Highest Tenth',
        'standard_name': 'sea_surface_wave_mean_period_of_highest_tenth',
        'units': 's',
        'comment': ('Wave period is the interval of time between repeated features on the waveform such as crests, '
                    'troughs or upward passes through the mean level. The mean wave period of highest tenth is the '
                    'mean period of the highest one-tenth of waves during the observation duration.'),
        'data_product_identifier': 'WAVSTAT-T10_L0',        
    },
    'wave_period_tp5': {
        'long_name': 'Peak Wave Period - Read Method',
        #'standard_name': 'sea_surface_wave_period_at_variance_spectral_density_maximum_read_method',
        'units': 's',
        'comment': ('Peak wave period in seconds as computed by the READ method. Tp5 has less statistical '
                    'variability than Tp because it is based on spectral moments. The Tp5 is determined '
                    'from calculating Fp5 which is the average frequency computed with the weighting '
                    'function S(f)**5 over the defined upper and lower frequency range.'),
        'data_product_identifier': 'WAVSTAT-TP5_L0',
    },
}


# +

def ratio_tp_to_hm0(tp, hm0):
    """
    This test evaluates the ratio of sig wave height to period.
    
    These tests are only applicable to values derived from the
    spectral analysis. Data products produced from the zero-
    crossing methods are unaffected. The results of this test
    are attributable to directional wave statistics products, i.e.
    those derived from the spectral variance values.
    
    :param tp: Peak wave period from the WAVSS (Tp) 
    :param hm0: Significant wave height calculated from spectral variance (hm0)
    :return qc_flags: a numpy array of flags indicating pass/fail (1/4) of the
        ratio test
    """
    # Find the breaks in the wave regime
    u10 = np.where(tp < 10)[0]
    o10 = np.where(tp >= 10)[0]
    
    # Next, we want to run the comparison on the low regime
    u10_good = (hm0[u10] >= 0.1)
    o10_good = (hm0[o10] >= 0.001*tp[o10]**2)

    # QC Check
    qc_flag = np.ones(tp.size)
    qc_flag[u10[~u10_good]] = 4
    qc_flag[o10[~o10_good]] = 4
    
    return qc_flag


def hsig_to_tavg(Hsig, Tavg):
    """
    This test evaluates wave heights against average wave periods.
    
    This test is based on empiriical fitting by NDBC. It tests wave
    height (Hsig) as a function of the average wave period to identify
    outliers. It is split into two regimes: one for Tavg < 5 sec and one
    for Tavg >=5 sec:
        * If Tavg < 5 sec: Hsig < 2.55 + (Tavg / 4)
        * If Tavg >= 5 sec: Hsig < (1.16 * Tavg) -2
    When Hsig exceeds those thresholds, it is flagged as suspicious.
    
    :param Hsig: Significant wave height calculated from the zero-down-crossing
        method
    :param Tavg: Average wave period
    :return qc_flag: a numpy array of flags indicating pass/suspicious (1/3)
    """
    # Find the breaks in the wave regime
    u5 = np.where(Tavg < 5)[0]
    o5 = np.where(Tavg >= 5)[0]
    
    # Run the comparison
    u5_good = (Hsig[u5] < (2.55 + (Tavg[u5] / 4)) )
    o5_good = (Hsig[o5] < ((1.16 * Tavg) - 2) )
    
    # Make the qc flags
    qc_flag = np.ones(Hsig.size)
    qc_flag[u5[~u5_good]] = 3
    qc_flag[o5[~o5_good]] = 3
    
    return qc_flag


def quality_checks(ds):
    """
    Assignment of QARTOD style quality flags to the WAVSS bulk wave 
    statistics. The two tests are the Tp-to-hm0 ratio test, which checks
    the data quality for parameters derived from spectral variance, and
    the Hsig-to-Tavg test, which test the significant wave height values.
    The assigned flag values are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail
        9 = Missing

    The final flag value represents the worst case assessment of the data
    quality.

    :param ds: xarray dataset with the WAVSS bulk wave statistics
    :return ds: dataset with QARTOD style quality flags added to the record
        per variable tested.
    """
    # Test the directional parameters quality
    directional_parameters = ["mean_direction", "mean_spectral_period", "mean_spread",
                              "peak_wave_period", "wave_height_hmo", "wave_period_tp5"]
    
    # Compute the directional parameter qc_flags
    directional_qc = ratio_tp_to_hm0(ds.peak_wave_period, ds.wave_height_hmo)
    
    # Add the directional quality_flags
    # add the qc_flags to the dataset, rolling up the results into a single value
    for p in directional_parameters:
        qc_summary = p + '_qc_summary_flag'
        if qc_summary in ds.variables:
            # add the new test results to the existing QC summary results
            qc = ds[qc_summary]
            flags = np.array([directional_qc, qc.values])
            ds[qc_summary] = ('time', flags.max(axis=0))
        else:
            # create a new QC summary variable
            ds[qc_summary] = ('time', directional_qc)

            # set up the attributes for the new variable
            ds[qc_summary].attrs = dict({
                'long_name': '%s QC Summary Flag' % ds[p].attrs['long_name'],
                'comment': ('Converts the QC Results values from a bitmap to a QARTOD style summary flag, where ',
                            'the values are 1 == pass, 2 == not evaluated, 3 == suspect or of high interest, ',
                            '4 == fail, and 9 == missing. The QC tests, as applied by OOI, only yield pass or ',
                            'fail values. By resetting, the QC flags become more user friendly and more nuanced.'),
                'flag_values': np.array([1, 2, 3, 4, 9]),
                'flag_meanings': 'pass not_evaluated suspect_or_of_high_interest fail missing'
            })

            # add the standard name if the variable has one
            if 'standard_name' in ds[p].attrs:
                ds[qc_summary].attrs['standard_name'] = '%s qc_summary_flag' % ds[p].attrs['standard_name']

    # Next, test the significant wave height
    wave_height_qc = hsig_to_tavg(ds.significant_wave_height, ds.peak_wave_period)
    qc_summary = "significant_wave_height" + '_qc_summary_flag'
    if qc_summary in ds.variables:
        # add the new test results to the existing QC summary results
        qc = ds[qc_summary]
        flags = np.array([wave_height_qc, qc.values])
        ds[qc_summary] = ('time', flags.max(axis=0))
    else:
        # create a new QC summary variable
        ds[qc_summary] = ('time', wave_height_qc)

        # set up the attributes for the new variable
        ds[qc_summary].attrs = dict({
            'long_name': '%s QC Summary Flag' % ds[p].attrs['long_name'],
            'comment': ('Converts the QC Results values from a bitmap to a QARTOD style summary flag, where ',
                        'the values are 1 == pass, 2 == not evaluated, 3 == suspect or of high interest, ',
                        '4 == fail, and 9 == missing. The QC tests, as applied by OOI, only yield pass or ',
                        'fail values. By resetting, the QC flags become more user friendly and more nuanced.'),
            'flag_values': np.array([1, 2, 3, 4, 9]),
            'flag_meanings': 'pass not_evaluated suspect_or_of_high_interest fail missing'
        })

        # add the standard name if the variable has one
        if 'standard_name' in ds[p].attrs:
            ds[qc_summary].attrs['standard_name'] = '%s qc_summary_flag' % ds[p].attrs['standard_name']
            
    return ds


# -

def wavss_datalogger(ds):
    """
    Takes WAVSS bulk wave statistics data recorded by the data loggers used in
    the CGSN/EA moorings and cleans up the data set to make it more user-friendly. 
    Primary task is renaming the alphabet soup parameter names and dropping some 
    parameters that are of no use/value.

    :param ds: initial WAVSS bulk wave statistics data set downloaded from OOI via the M2M system
    :return ds: cleaned up data set
    """
    drop_list = ['dcl_controller_timestamp', 'internal_timestamp', 'time_string',
                 'date_string', 'serial_number']
    for var in ds.variables:
        if var in drop_list:
            ds = ds.drop(var)

    # reset some attributes
    for key, value in ATTRS.items():
        for atk, atv in value.items():
            if key in ds.variables:
                ds[key].attrs[atk] = atv

    # parse the OOI QC variables and add QARTOD style QC summary flags to the data, converting the
    # bitmap represented flags into an integer value representing pass == 1, suspect or of high
    # interest == 3, and fail == 4.
    ds = parse_qc(ds)
    
    # Run and add in the additional quality checks
    ds = quality_checks(ds)

    return ds
