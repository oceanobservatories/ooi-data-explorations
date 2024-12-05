import numpy as np
import pandas as pd
import xarray as xr

QCThresholds = {
    'error_velocity': {
        'pass': 0.05,
        'fail': 0.20 },
    'correlation_magnitude': {
        'pass': 115,
        'fail': 63 },
    'vertical_velocity': {
        'pass': 0.30,
        'fail': 0.50 },
    'horizontal_speed': {
        'pass': 1.00,
        'fail': 2.50 },
    'percent_good': {
        'ADCPT': {
            'pass': 56,
            'fail': 45 },
        'ADCPS': {
            'pass': 48,
            'fail': 38 }
    }
}


def sidelobe_depth(ds, theta=20):
    """
    Calculate the sidelobe contamination depth for the given ADCP.

    The sidelobe intereference depth 
    is caluclated following Lentz et al. (2022) where:
    
        z_ic = ha*[1 - cos(theta)] + 3*delta_Z/2
        
    z_ic is the depth above which there is sidelobe interference,
    ha is the transducer face depth, theta is the beam angle, and
    delta_Z is the cell-bin depth. At this time we ignore instrument 
    tilt at this time and its impact, assuming fixed beam angle.

    Parameters
    ----------
    ds: 
    theta: transducer head angle (default 20 degrees)

    Returns
    -------
    z_ic: the maximum depth impacted by sidelobe interference (m)
    """
    # First, get the transducer depth
    depth = ds['depth_from_pressure']
    ha = depth.interpolate_na(dim='time', method='linear')

    # Next, get the beam angle
    theta = np.deg2rad(theta)

    # Grab the cell length and convert to m
    delta_z = ds['cell_length'].mean(skipna=True)/100

    # Calculate the range of cells contaminated by sidelobe interference
    z_ic = ha*(1 - np.cos(theta)) + 3*delta_z/2

    return z_ic


def sidelobe_qc(ds):
    """
    Add sidelobe interference QARTOD-style flags to the ADCP dataset

    Assignment of QARTOD style quality flags to the ADCP velocity data based
    on estimation of sidelobe contamination. The sidelobe intereference depth 
    is caluclated following Lentz et al. (2022) where:
    
        z_ic = ha*[1 - cos(theta)] + 3*delta_Z/2
        
    z_ic is the depth above which there is sidelobe interference,
    ha is the transducer face depth, theta is the beam angle, and
    delta_Z is the cell-bin depth. Bin depths less than z_ic are
    considered contaminated and the data is considered bad. The assigned flag
    values are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail
        9 = Missing

    Parameters
    ----------
    ds: xarray dataset of downloaded adcp data from OOINet

    Returns
    -------
    ds: xarray dataset with QC parameter for sidelobe interference
    """
    # First, get the sidelobe contamination depth z_ic
    z_ic = sidelobe_depth(ds)

    # Next, create a qc_flag for each bin measurement and
    # identify the bins which are shallower than z_ic
    qc_flag = np.ones(ds['bin_depths'].shape, dtype=int)
    mask = ds['bin_depths'] < z_ic
    qc_flag[mask] = 4

    # Add the qc_flags
    qc_name = 'bin_depths_qc_summary_flag'
    ds[qc_name] =  (['time', 'bins'], qc_flag)
    
     # set up the attributes for the new variable
    ds[qc_name].attrs = dict({
        'long_name': '%s QC Summary Flag' % ds['bin_depths'].attrs['long_name'],
        'comment': ('A QARTOD style summary flag indicating depth bins with sidelobe contamination, where ',
                    'the values are 1 == pass, 2 == not evaluated, 3 == suspect or of high interest, ',
                    '4 == fail, and 9 == missing. The QC tests, as applied by OOI, only yield pass or ',
                    'fail values. Sidelobe contamination depth is determined following Lentz et al (2022).'),
        'flag_values': np.array([1, 2, 3, 4, 9]),
        'flag_meanings': 'pass not_evaluated suspect_or_of_high_interest fail missing'
    })

    return ds


def error_velocity_qc(ds, suspect, fail):
    """
    Determine ADCP QC based on Error velocity. This algorithm
    uses thresholds from the NDBC QC Handbook based on the 
    kHz of the relevant TRDI ADCP.
    """
    # Set up a qc_flags the shape of the variable
    qc_flags = np.ones(ds['error_seawater_velocity'].shape, dtype=int)
    
    # Now find the "suspect" values
    mask = (np.abs(ds['error_seawater_velocity']) > suspect)
    qc_flags[mask] = 3
    
    # Now find the "fail" values
    mask = (np.abs(ds['error_seawater_velocity']) > fail)
    qc_flags[mask] = 4

    return qc_flags


def vertical_velocity_qc(ds, suspect, fail):
    """
    Determine ADCP QC based on vertical velocity. This algorithm
    uses thresholds from the NDBC QC Handbook based on the 
    kHz of the relevant TRDI ADCP.
    """
    # Set up a qc_flags the shape of the variable
    qc_flags = np.ones(ds['upward_seawater_velocity'].shape, dtype=int)
    
    # Now find the "suspect" values
    mask = (np.abs(ds['upward_seawater_velocity']) > suspect)
    qc_flags[mask] = 3
    
    # Now find the "fail" values
    mask = (np.abs(ds['upward_seawater_velocity']) > fail)
    qc_flags[mask] = 4

    return qc_flags


def horizontal_speed_qc(ds, suspect, fail):
    """
    Determine ADCP QC based on vertical velocity. This algorithm
    uses thresholds from the NDBC QC Handbook based on the 
    kHz of the relevant TRDI ADCP.

    The good, suspect, fail are defined as follows:
        good: both east and north velocities are good OR one is suspect
              and the other is good
        suspect: both east and north velocities are suspect
        fail: either east or north velocities fail
    """
    # Set up a qc_flags the shape of the variable
    qc_flags_east = np.ones(ds['eastward_seawater_velocity'].shape, dtype=int)
    qc_flags_north = np.ones(ds['northward_seawater_velocity'].shape, dtype=int)
    
    # Now find the "suspect" values
    mask = (np.abs(ds['eastward_seawater_velocity']) > suspect)
    qc_flags_east[mask] = 3
    mask = (np.abs(ds['northward_seawater_velocity']) > suspect)
    qc_flags_north[mask] = 3
    
    # Now find the "fail" values
    mask = (np.abs(ds['eastward_seawater_velocity']) > fail)
    qc_flags_east[mask] = 4
    mask = (np.abs(ds['northward_seawater_velocity']) > fail)
    qc_flags_north[mask] = 4

    # Now combine them using math to parse out the combinations
    qc_flags = np.ones(ds['eastward_seawater_velocity'].shape, dtype=int)
    
    # Suspect flags when both directions are suspect
    suspect_flags = ((qc_flags_east == 3) & (qc_flags_north == 3))
    qc_flags[suspect_flags] = 3

    # Fail flags when either direction is bad
    bad_flags = ((qc_flags_east == 4) | (qc_flags_north == 4))
    qc_flags[bad_flags] = 4

    return qc_flags


def correlation_magnitude_qc(ds, suspect, fail):
    """
    Determine ADCP QC based on correlation magnitude. This algorithm
    uses thresholds from the NDBC QC Handbook based on the 
    kHz of the relevant TRDI ADCP.

    The good, suspect, fail are defined as follows:
        good: correlation magnitudes of at least 3 out of 4 beams pass
        suspect: only two of the beams pass
        fail: one or none of the beams pass
    """
    # Set up the qc_flags 
    qc_flags = np.ones(ds['correlation_magnitude_beam1'].shape, dtype=int)
    
    # Simplify implementation by adding booleans
    beam1_pass = (ds['correlation_magnitude_beam1'] > suspect).astype(int)
    beam2_pass = (ds['correlation_magnitude_beam1'] > suspect).astype(int)
    beam3_pass = (ds['correlation_magnitude_beam1'] > suspect).astype(int)
    beam4_pass = (ds['correlation_magnitude_beam1'] > suspect).astype(int)

    # Sum the results
    total_pass = beam1_pass + beam2_pass + beam3_pass + beam4_pass

    # Good values sum to 3 or 4
    # Suspect values sum to 2
    mask = (total_pass == 2)
    qc_flags[mask] = 3

    # Fail values sum to 0 or 1
    mask = (total_pass < 2)
    qc_flags[mask] = 4

    # Return the qc_flags
    return qc_flags


def percent_good_qc(ds, suspect, fail):
    """
    Determine ADCP QC based on correlation magnitude. This algorithm
    uses thresholds from the NDBC QC Handbook based on the 
    kHz of the relevant TRDI ADCP. Percent good is calculated from the
    best returns of the 3 and 4 beam solutions, since either may be
    used to calculate the seawater velocities.
    """
    # Merge the 3 and 4 beam solutions along a new axis and take the maximum
    percent_good = np.stack([ds['percent_good_3beam'].values, ds['percent_good_4beam'].values], axis=-1)
    percent_good = np.max(percent_good, axis=2)
    
    # Create a qc_flags array
    qc_flags = np.ones(percent_good.shape, dtype=int)
    
    # Find where the flags are suspect
    mask = (percent_good < suspect)
    qc_flags[mask] = 3

    # Find where the correlation magnitude is bad
    mask = (percent_good < fail)
    qc_flags[mask] = 4

    # Return the results
    return qc_flags


def merge_qc(test_results):
    """
    Merge the results of the different QC tests into a single
    output. The results for the entire ensemble are:
        Pass:    100% of qc tests pass
        Suspect: At least 50% of tests pass or are suspect
        Fail:    Less than 50% of tests pass or are suspect
    """
    n = len(test_results)
    qc_flags = np.zeros(test_results[0].shape, dtype=int)

    # First, calculate the most "inclusive" case of the number
    # of pass OR suspect tests and calculate the fraction 
    tests_suspect = np.zeros(test_results[0].shape, dtype=int)
    for test in test_results:
        suspect = ((test == 1) | (test == 3)).astype(int)
        tests_suspect = tests_suspect + suspect
    # Now calculate the fraction of tests that passed
    total_suspect = (tests_suspect / n)
    # Find where not enough tests pass and mark as fail
    mask = (total_suspect < 0.5)
    qc_flags[mask] = 4
    # Mark the rest as suspect. Will test for pass next
    mask = (total_suspect >= 0.5)
    qc_flags[mask] = 3

    # Now test for if all tests pass
    tests_passed = np.zeros(test_results[0].shape, dtype=int)
    # Sum up the number of "passes" in the qc_tests
    for test in test_results:
        passed = (test == 1).astype(int)
        tests_passed = tests_passed + passed
    # Now calculate the fraction of tests that passed
    total_passed = (tests_passed / n)
    mask = (total_passed == 1)
    qc_flags[mask] = 1

    # Return the result
    return qc_flags
