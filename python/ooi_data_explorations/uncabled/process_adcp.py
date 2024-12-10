import numpy as np
import numpy.typing as npt
import pandas as pd
import xarray as xr

# Can use these based on NDBC but recommend using the 
# TRDI QAQC Model rev12-1
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


def sidelobe_depth(ds: xr.Dataset, theta: int = 20) -> xr.DataArray:
    """
    Calculate the sidelobe contamination depth for the given ADCP.

    The sidelobe intereference depth 
    is caluclated following Lentz et al. (2022) where:
    
        z_ic = ha*[1 - cos(theta)] + 3*delta_Z/2
        
    z_ic is the depth above which there is sidelobe interference,
    ha is the transducer face depth, theta is the beam angle, and
    delta_Z is the cell-bin depth. We ignore instrument tilt at 
    this time and its impact, assuming fixed beam angle.

    Parameters
    ----------
    ds: xarray.dataset
        A TRDI ADCP dataset from OOI from which to calculate the
        sidelobe interference depth
    theta: float, Default = 20
        Beam angle of the given ADCP

    Returns
    -------
    z_ic:
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


def sidelobe_qc(ds: xr.Dataset) -> xr.Dataset:
    """
    Add sidelobe interference QARTOD-style flags to the ADCP dataset

    Assignment of QARTOD style quality flags to the ADCP velocity data based
    on estimation of sidelobe contamination. The sidelobe intereference depth 
    is calculated following Lentz et al. (2022) where:
    
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
    ds: xarray.Dataset
        The dataset containing the TRDI ADCP data downloaded from
        OOINet in .netcdf format

    Returns
    -------
    ds: xarray.Dataset
        The input dataset with the sidelobe interference QC flag
        added as a dataset variable
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


def error_velocity_qc(ds: xr.Dataset, suspect: float | int, fail: float | int) -> npt.NDArray[int]:
    """
    Determine ADCP QC based on Error velocity and assign
    QARTOD-style flags. This algorithm uses thresholds computed
    using the TRDI ADCP Data AQ-QC Model rev12-1. The assigned 
    flag values are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail

    The pass, suspect, fail flags are defined as follows:
        pass: error velocity is less than the suspect threshold
        suspect: error velocity is between the suspect and fail thresholds 
        fail: error velocity exceeds the fail threshold

    Parameters
    ----------
    ds: xarray.Dataset
        The dataset containing the TRDI ADCP data downloaded from
        OOINet in .netcdf format
    suspect: float
        The suspect threshold computed from the TRDI QA-QC Model
    fail: flot
        The fail threshold computed from the TRDI QA-QC Model

    Returns
    -------
    qc_flags: numpy.array[int]
        An array of QARTOD-style flags indicating the results of the QC
        test for each given datum

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


def vertical_velocity_qc(ds: xr.Dataset, suspect: float | int, fail: float | int) -> npt.NDArray[int]:
    """
    Determine ADCP QC based on vertical velocity and assign
    QARTOD-style flags. This algorithm uses thresholds computed
    using the TRDI ADCP Data AQ-QC Model rev12-1. The assigned 
    flag values are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail

    The pass, suspect, fail flags are defined as follows:
        pass: vertical velocity is less than the suspect threshold
        suspect: vertical velocity is between the suspect and fail thresholds 
        fail: vertical velocity exceeds the fail threshold

    Parameters
    ----------
    ds: xarray.Dataset
        The dataset containing the TRDI ADCP data downloaded from
        OOINet in .netcdf format
    suspect: float
        The suspect threshold computed from the TRDI QA-QC Model
    fail: flot
        The fail threshold computed from the TRDI QA-QC Model

    Returns
    -------
    qc_flags: numpy.array[int]
        An array of QARTOD-style flags indicating the results of the QC
        test for each given datum
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


def horizontal_speed_qc(ds: xr.Dataset, suspect: float, fail: float) -> npt.NDArray[int]:
    """
    Determine ADCP QC based on vertical velocity and assign
    QARTOD-style flags. This algorithm uses thresholds computed
    using the TRDI ADCP Data AQ-QC Model rev12-1. The assigned 
    flag values are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail

    The pass, suspect, fail flags are defined as follows:
        pass: both east and north velocities are good OR one is suspect
              and the other is good
        suspect: both east and north velocities are suspect
        fail: either east or north velocities fail

    Parameters
    ----------
    ds: xarray.Dataset
        The dataset containing the TRDI ADCP data downloaded from
        OOINet in .netcdf format
    suspect: float
        The suspect threshold computed from the TRDI QA-QC Model
    fail: flot
        The fail threshold computed from the TRDI QA-QC Model

    Returns
    -------
    qc_flags: numpy.array[int]
        An array of QARTOD-style flags indicating the results of the QC
        test for each given datum
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


def correlation_magnitude_qc(ds: xr.Dataset, suspect: float, fail: float) -> npt.NDArray[int]:
    """
    Determine ADCP QC based on correlation magnitude and assign
    QARTOD-style flags. This algorithm uses thresholds computed
    using the TRDI ADCP Data AQ-QC Model rev12-1. The assigned 
    flag values are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail

    The pass, suspect, fail are defined as follows:
        pass: correlation magnitudes of at least 3 out of 4 beams pass
        suspect: only two of the beams pass
        fail: one or none of the beams pass

    Parameters
    ----------
    ds: xarray.Dataset
        The dataset containing the TRDI ADCP data downloaded from
        OOINet in .netcdf format
    suspect: float
        The suspect threshold computed from the TRDI QA-QC Model
    fail: flot
        The fail threshold computed from the TRDI QA-QC Model

    Returns
    -------
    qc_flags: numpy.array[int]
        An array of QARTOD-style flags indicating the results of the QC
        test for each given datum
    """
    # Set up the qc_flags 
    qc_flags = np.ones(ds['correlation_magnitude_beam1'].shape, dtype=int)
    
    # Simplify implementation by adding booleans
    beam1_pass = (ds['correlation_magnitude_beam1'] > suspect).astype(int)
    beam2_pass = (ds['correlation_magnitude_beam2'] > suspect).astype(int)
    beam3_pass = (ds['correlation_magnitude_beam3'] > suspect).astype(int)
    beam4_pass = (ds['correlation_magnitude_beam4'] > suspect).astype(int)

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


def percent_good_qc(ds: xarray.Dataset, suspect: float, fail: float) -> npt.NDArray[int]:
    """
    Determine ADCP QC based on the percent good returned for each
    beam and assign QARTOD-style flags. This algorithm uses thresholds
    computed using the TRDI ADCP Data AQ-QC Model rev12-1. The assigned 
    flag values are:

        1 = Pass
        3 = Suspect or of High Interest
        4 = Fail. This algorithm

    Percent good is calculated from the best returns of either 3 and 4 beam
    solutions, since either may be used to calculate the seawater velocities.

    Parameters
    ----------
    ds: xarray.Dataset
        The dataset containing the TRDI ADCP data downloaded from
        OOINet in .netcdf format
    suspect: float
        The suspect threshold computed from the TRDI QA-QC Model
    fail: flot
        The fail threshold computed from the TRDI QA-QC Model

    Returns
    -------
    qc_flags: numpy.array[int]
        An array of QARTOD-style flags indicating the results of the QC
        test for each given datum
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


def merge_qc(test_results: list[npt.NDArray[int]]) -> npt.NDArray[int]:
    """
    Merge the results of the different QC tests into a single
    output. The results for the entire ensemble are:
    
        Pass:    100% of qc tests pass
        Suspect: At least 50% of tests pass or are suspect
        Fail:    Less than 50% of tests pass or are suspect

    The assigned QARTOD-style flag values are:
        1 = Pass
        3 = Suspect
        4 = Fail 

    Parameters
    ----------
    test_results: list[numpy.array[int]]
        A list containing all of the run individual ADCP
        QC test results

    Returns
    -------
    qc_flags: numpy.array[int]
        A numpy array that contains the combined results
        of the individual QC tests passed in test_results
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
