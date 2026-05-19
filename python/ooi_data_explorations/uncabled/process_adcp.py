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


def _threshold_qc(values: npt.NDArray, suspect: float, fail: float) -> npt.NDArray[np.int8]:
    """
    Apply simple absolute-value threshold QC and return QARTOD flags.

    Uses np.select for a single, vectorised pass over the data rather than
    applying two sequential boolean masks.

    Parameters
    ----------
    values: numpy.ndarray
        Raw measurement array (already extracted from the dataset).
    suspect: float
        Threshold above which a value is considered suspect.
    fail: float
        Threshold above which a value is considered a fail.

    Returns
    -------
    qc_flags: numpy.ndarray[int8]
    """
    abs_values = np.abs(values)
    return np.select(
        [abs_values > fail, abs_values > suspect],
        [4, 3],
        default=1,
    ).astype(np.int8)


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
    ha = ds['depth_from_pressure'].interpolate_na(dim='time', method='linear')

    # Next, get the beam angle
    theta = np.deg2rad(theta)

    # Grab the cell length and convert from cm to m
    delta_z = ds['cell_length'].mean(skipna=True) / 100

    # Calculate the range of cells contaminated by sidelobe interference
    z_ic = ha * (1 - np.cos(theta)) + 3 * delta_z / 2
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

    # Identify bins which are shallower than z_ic
    # QC flag: start with all-pass, mark contaminated bins as fail
    qc_flag = np.where(ds['bin_depths'] < z_ic, 4, 1).astype(np.int8)

    # Add the qc_flags and associated attributes
    qc_name = 'bin_depths_qc_summary_flag'
    ds[qc_name] =  (['time', 'bins'], qc_flag)    
    ds[qc_name].attrs = {
        'long_name': '%s QC Summary Flag' % ds['bin_depths'].attrs['long_name'],
        'comment': ('A QARTOD style summary flag indicating depth bins with sidelobe contamination, where ',
                    'the values are 1 == pass, 2 == not evaluated, 3 == suspect or of high interest, ',
                    '4 == fail, and 9 == missing. The QC tests, as applied by OOI, only yield pass or ',
                    'fail values. Sidelobe contamination depth is determined following Lentz et al (2022).'),
        'flag_values': np.array([1, 2, 3, 4, 9]),
        'flag_meanings': 'pass not_evaluated suspect_or_of_high_interest fail missing'
    }
    return ds


def error_velocity_qc(ds: xr.Dataset, suspect: float | int, fail: float | int) -> npt.NDArray[np.int8]:
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
    return _threshold_qc(ds['error_seawater_velocity'].values, suspect, fail)


def vertical_velocity_qc(ds: xr.Dataset, suspect: float | int, fail: float | int) -> npt.NDArray[np.int8]:
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
    return _threshold_qc(ds['upward_seawater_velocity'].values, suspect, fail)


def horizontal_speed_qc(ds: xr.Dataset, suspect: float, fail: float) -> npt.NDArray[np.int8]:
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
    # Get the east and north velocity
    east = np.abs(ds['eastward_seawater_velocity'].values)
    north = np.abs(ds['northward_seawater_velocity'].values)
    
    # Fail if either component exceeds the fail threshold; suspect if both
    # exceed the suspect threshold; otherwise pass — single vectorised call.
    qc_flags = np.select(
        [(east > fail) | (north > fail),
         (east > suspect) & (north > suspect)],
        [4, 3],
        default=1,
    ).astype(np.int8)

    return qc_flags


def correlation_magnitude_qc(ds: xr.Dataset, suspect: float, fail: float) -> npt.NDArray[np.int8]:
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
    # Stack all four beams along a new trailing axis, then count passes in one op
    beams = np.stack([
        ds['correlation_magnitude_beam1'].values,
        ds['correlation_magnitude_beam2'].values,
        ds['correlation_magnitude_beam3'].values,
        ds['correlation_magnitude_beam4'].values,
    ], axis=-1)

    # Get the total pass
    total_pass = np.sum(beams > suspect, axis=-1)  # shape: (time, bins)

    # Find which flags passes
    qc_flags = np.select(
        [total_pass < 2, total_pass == 2],
        [4, 3],
        default=1,
    ).astype(np.int8)

    # Return the qc_flags
    return qc_flags


def percent_good_qc(ds: xr.Dataset, suspect: float, fail: float) -> npt.NDArray[np.int8]:
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
    # Take the element-wise maximum of the 3-beam and 4-beam solutions
    percent_good = np.maximum(
        ds['percent_good_3beam'].values,
        ds['percent_good_4beam'].values,
    )

    # Get the qc_flags
    qc_flags = np.select(
        [percent_good < fail, percent_good < suspect],
        [4, 3],
        default=1,
    ).astype(np.int8)

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
    # Stack into (n_tests, ...) then compute per-element statistics in one pass
    stacked = np.stack(test_results, axis=0)  # shape: (n, time, bins)
    n = len(test_results)

    # Count passes and pass-or-suspects across all tests simultaneously
    n_passed   = np.sum(stacked == 1, axis=0)
    n_suspect  = np.sum((stacked == 1) | (stacked == 3), axis=0)

    # Calculate the qc_flags
    qc_flags = np.select(
        [n_suspect / n < 0.5,   # fail: fewer than half are even suspect
         n_passed  / n < 1.0],  # suspect: not all passed
        [4, 3],
        default=1,              # pass: every test passed
    ).astype(np.int8)

    # Return the result
    return qc_flags