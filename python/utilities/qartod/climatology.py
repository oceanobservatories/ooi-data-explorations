import numpy as np
import pandas as pd
import xarray as xr


def calc_regression_climatology(time_series, freq=1/12, lin_trend=False):
    """
    Calculate a two-cycle harmonic linear regression following Ax=b.

    This is an Ordinary-Least Squares regression fit for a two-cycle harmonic
    for OOI climatology data.

    Parameters
    ----------
    time_series: (numpy.array)
        A numpy array of monthly-binned mean data to be fitted with a harmonic
            cycle.
    freq: (float)
        The frequency of the fit (default = 1/12).
    lin_trend: (boolean)
        Switch to determine if a linear trends should be added to the
        climatological fit (default=False).

    Returns
    -------
    seasonal_cycle: (numpy.array)
        A numpy array of the monthly-best fit values fitted with OLS-regressed
            harmonic cycle. Note this is NOT robust to significant outliers.
    beta: (numpy.array)
        A numpy array of the coefficients of best-fit
    sigma: (float)
        The standard deviation of the monthly-best fit values against the input
            time series.
    """
    # Rename some of the data variables
    ts = time_series
    N = len(ts)
    t = np.arange(0, N, 1)
    new_t = t
    f = freq

    # Drop NaNs from the fit
    mask = np.isnan(ts)
    ts = ts[mask == False]
    t = t[mask == False]
    N = len(t)

    # Build the linear coefficients as a stacked array (this is the matrix A)
    if lin_trend:
        arr0 = np.ones(N)
        arr1 = np.sin(2*np.pi*f*t)
        arr2 = np.cos(2*np.pi*f*t)
        arr3 = np.sin(4*np.pi*f*t)
        arr4 = np.cos(4*np.pi*f*t)
        x = np.stack([arr0, arr1, arr2, arr3, arr4, t])
    else:
        arr0 = np.ones(N)
        arr1 = np.sin(2*np.pi*f*t)
        arr2 = np.cos(2*np.pi*f*t)
        arr3 = np.sin(4*np.pi*f*t)
        arr4 = np.cos(4*np.pi*f*t)
        x = np.stack([arr0, arr1, arr2, arr3, arr4])

    # Fit the coefficients using OLS
    beta, _, _, _ = np.linalg.lstsq(x.T, ts)

    # Now fit a new timeseries with the coefficients of best fit
    if lin_trend:
        seasonal_cycle = beta[0] + beta[1]*np.sin(2*np.pi*f*new_t)
        + beta[2]*np.cos(2*np.pi*f*new_t) + beta[3]*np.sin(4*np.pi*f*new_t)
        + beta[4]*np.cos(4*np.pi*f*new_t) + beta[-1]*new_t
    else:
        seasonal_cycle = beta[0] + beta[1]*np.sin(2*np.pi*f*new_t)
        + beta[2]*np.cos(2*np.pi*f*new_t) + beta[3]*np.sin(4*np.pi*f*new_t)
        + beta[4]*np.cos(4*np.pi*f*new_t)

    # Now calculate the standard deviation of the time series
    sigma = np.sqrt((1/(len(ts)-1))
                    * np.sum(np.square(ts - seasonal_cycle[mask == False])))

    return seasonal_cycle, beta, sigma


def qartod_climatology(ds):
    """
    Calculate the monthly QARTOD Climatology for a time series.

    Parameters
    ----------
    ds: (xarray.dataArray)
        An xarray data array with main dimension of time.

    Returns
    -------
    results: (list of tuples)
        A list of tuples in the format of
        (numerical month, None, [lower bound, upper bound], month)
    """
    # Calculate the monthly means of the dataset
    monthly = ds.resample(time="M").mean()

    # Fit the regression for the monthly harmonic
    cycle, beta, sigma = calc_regression_climatology(monthly.values)

    # Calculate the monthly means, take a look at the seasonal cycle values
    climatology = pd.Series(cycle, index=monthly.time.values)
    climatology = climatology.groupby(climatology.index.month).mean()

    # Now add the standard deviations to get the range of data
    lower = np.round(climatology-sigma*2, decimals=2)
    upper = np.round(climatology+sigma*2, decimals=2)

    # This generates the results tuple
    results = []
    for month in climatology.index:
        tup = (month, None, [lower[month], upper[month]], None)
        results.append(tup)

    return results
