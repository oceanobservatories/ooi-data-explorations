import os
import numpy as np
import pandas as pd
import xarray as xr
import re


class Climatology():
    """Climatology fitting process for QARTOD.

    For a given parameter in a dataset, fit a climatology cycle using either
    montly or daily binning to the data, and format the data into a qcConfig
    object optimized for use with Axiom-implemented QARTOD climatology test.

    Example
    -------
    from qartod.climatology import Climatology
    climatology = Climatology()
    climatology.fit(pco2_dataset, "pco2_seawater")
    climatology.make_qcConfig()
    climatology.qcConfig = {'qartod': {'climatology': {'config':
        [{'tspan': [0, 1], 'vspan': [347.24, 511.88], 'period': 'month'},
         {'tspan': [1, 2], 'vspan': [321.59, 486.23], 'period': 'month'},
         {'tspan': [2, 3], 'vspan': [313.65, 478.29], 'period': 'month'},
         {'tspan': [3, 4], 'vspan': [325.54, 490.18], 'period': 'month'},
         {'tspan': [4, 5], 'vspan': [354.08, 518.72], 'period': 'month'},
         {'tspan': [5, 6], 'vspan': [391.63, 556.27], 'period': 'month'},
         {'tspan': [6, 7], 'vspan': [428.12, 592.76], 'period': 'month'},
         {'tspan': [7, 8], 'vspan': [453.77, 618.41], 'period': 'month'},
         {'tspan': [8, 9], 'vspan': [461.71, 626.35], 'period': 'month'},
         {'tspan': [9, 10], 'vspan': [449.82, 614.46], 'period': 'month'},
         {'tspan': [10, 11], 'vspan': [421.28, 585.92], 'period': 'month'},
         {'tspan': [11, 12], 'vspan': [383.73, 548.37], 'period': 'month'}]}}}"""

    def resample(self, ds, param, period):
        """Resample a data variable from the dataset to the desired period.

        Parameters
        ----------
        ds: (xarray.DataSet)
            An xarray dataset containing the given data variable to be
            resampled and a primary dimenstion of time.
        param: (string)
            The name of the data variable to resample from the given dataset.
        period: (string)
            The time period (e.g. month = "M", day = "D") to resample and take
            the mean of.

        Returns
        -------
        da: (xarray.DataArray)
            An xarray DataArray with the resampled mean values for the given
            param.
        """
        df = ds[param].to_dataframe()
        da = xr.DataArray(df.resample(period).mean())

        return da

    def fit(self, ds, param, period="M", cycles=1, lin_trend=False):
        """Fit the climatology with either monthly or daily binning.

        Parameters
        ----------
        ds: (xarray.DataSet)
            An xarray datasets containing the given data variable to be fitted,
            with a primary dimension of time.
        param: (string)
            The name of the data variable from the given dataset to fit
        period: (string)
            The time period (e.g. month = "M", day = "D") to bin the data,
            which will correspond to the fitted result
        cycle: (int: 1)
            The number of cycles per year to fit the data with
        lin_trend: (bool: False)
            Whether to include a monotonic linear trend in the fitted data
        """
        # Resample the data
        da = self.resample(ds, param, period)

        # Calculate the frequency from the period
        if period == "M":
            freq = 1/12
        elif period == "D":
            freq = 1/365
        else:
            pass

        # Get the time series
        time_series = da.values.reshape(-1)

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

        arr0 = np.ones(N)
        if cycles == 1:
            arr1 = np.sin(2*np.pi*f*t)
            arr2 = np.cos(2*np.pi*f*t)
            if lin_trend:
                x = np.stack([arr0, arr1, arr2, t])
            else:
                x = np.stack([arr0, arr1, arr2])
        else:
            arr1 = np.sin(2*np.pi*f*t)
            arr2 = np.cos(2*np.pi*f*t)
            arr3 = np.sin(4*np.pi*f*t)
            arr4 = np.cos(4*np.pi*f*t)
            if lin_trend:
                x = np.stack([arr0, arr1, arr2, arr3, arr4, t])
            else:
                x = np.stack([arr0, arr1, arr2, arr3, arr4])

        # Fit the coefficients using OLS
        beta, _, _, _ = np.linalg.lstsq(x.T, ts)

        # Now fit a new timeseries with the coefficients of best fit
        if cycles == 1:
            if lin_trend:
                fitted_data = beta[0] + beta[1]*np.sin(2*np.pi*f*new_t)
                + beta[2]*np.cos(2*np.pi*f*new_t)
                + beta[-1]*new_t
            else:
                fitted_data = beta[0] + beta[1]*np.sin(2*np.pi*f*new_t)
                + beta[2]*np.cos(2*np.pi*f*new_t)
        else:
            if lin_trend:
                fitted_data = beta[0] + beta[1]*np.sin(2*np.pi*f*new_t)
                + beta[2]*np.cos(2*np.pi*f*new_t)
                + beta[3]*np.sin(4*np.pi*f*new_t)
                + beta[4]*np.cos(4*np.pi*f*new_t)
                + beta[-1]*new_t
            else:
                fitted_data = beta[0] + beta[1]*np.sin(2*np.pi*f*new_t)
                + beta[2]*np.cos(2*np.pi*f*new_t)
                + beta[3]*np.sin(4*np.pi*f*new_t)
                + beta[4]*np.cos(4*np.pi*f*new_t)

        # Now calculate the standard deviation of the time series
        sigma = np.sqrt((1/(len(ts)-1))*np.sum(np.square(ts - fitted_data[mask == False])))
        sigma = np.round(sigma, decimals=2)

        # Reformat the fitted data into a pandas series indexed by the time and
        # store the period information
        fitted_data = pd.Series(data=np.round(fitted_data, decimals=2),
                                index=da.time.values)
        fitted_data.index.freq = period

        # Reformat
        beta = np.round(beta, decimals=2)

        # Save the results as attributes of the object
        self.fitted_data = fitted_data
        self.sigma = sigma
        self.beta = beta

    def make_config(self):
        """Calculate the config dictionary for climatology."""
        config = []

        months = np.arange(1, 13, 1)

        for month in months:
            val = self.fitted_data[self.fitted_data.index.month == month]
            if len(val) == 0:
                val = np.nan
            else:
                val = val.mean()

            # Get the min/max values
            vmin = np.round(val-self.sigma*3, 2)
            vmax = np.round(val+self.sigma*3, 2)

            # Record the results
            tspan = [month-1, month]
            vspan = [vmin, vmax]

            # Add in the tspan, vspan, and period information
            config.append({
                "tspan": tspan,
                "vspan": vspan,
                "period": "month"
            })

        return config

    def make_qcConfig(self):
        """Build properly formatted qcConfig object for qartod climatology."""
        config = {
            "qartod": {
                "climatology": {
                    "config": self.make_config()
                }
            }
        }

        self.qcConfig = config
