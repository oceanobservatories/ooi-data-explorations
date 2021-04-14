import numpy as np
import pandas as pd
import xarray as xr


class Climatology():
    """Climatology fitting process for QARTOD.

    For a given parameter in a dataset, fit a climatology cycle using either
    montly or daily binning to the data, and format the data into a qcConfig
    object optimized for use with Axiom-implemented QARTOD climatology test.

    Example
    -------
    from qartod.climatology import Climatology
    climatology = Climatology()
    climatology.fit(ctdbp_dataset, "ctdbp_seawater_temperature")
    """

    def _std(self, x):
        """Calculate the grouped standard deviations."""
        N = len(x)
        std = np.sqrt(np.sum(x**2)/N)
        return std

    def calc_standard_deviations(self, ds, param):
        """Calculate the monthly standard deviations.

        Calculate the standard deviations from the original observations
        and the monthly expectation derived from the two-cycle harmonic
        fit to the monthly mean values.

        Parameters
        ----------
        ds: (xarray.DataSet)
            DataSet of the original time series observations
        param: (str)
            A string corresponding to the variable in the DataSet to fit.

        Attributes
        -------
        monthly_avg: (pandas.Series)
            The climatological expectation for each calendar month of a year
        monthly_std: (pandas.Series)
            The standard deviation calculated from the observations and
            the climatological expectation for each calendard month of a year
        """
        # Monthly expectation mu
        mu = xr.DataArray(self.fitted_data)
        mu = mu.groupby(mu.time.dt.month).mean()

        # Group the original observations by month
        X = ds[param].groupby(ds.time.dt.month)

        # Calculate the difference between the obs and expectation
        diff = X - mu

        # Calculate the standard deviations for each month
        std = diff.groupby("month").apply(self._std)

        # Save the results
        self.monthly_avg = pd.Series(mu.values, index=mu.month)
        self.monthly_std = pd.Series(std.values, index=std.month)

    def fit(self, ds, param):
        """Calculate the climatological fit and monthly standard deviations.

        Calculates the climatological fit for a time series. First, the data
        are binned by month and averaged. Next, a two-cycle harmonic is fitted
        via OLS-regression. The climatological expected value for each month
        is then calculated from the regression coefficients. Finally, the
        standard deviation is derived using the observations for a given month
        and the climatological fit for that month as the expected value.

        Parameters
        ----------
        ds: (xarray.DataSet)
            DataSet of the original time series observations
        param: (str)
            A string corresponding to the variable in the DataSet to fit.

        Attributes
        -------
        fitted_data: (pandas.Series)
            The climatological monthly expectation calculated from the
            regression, indexed by the year-month
        regression: (dict)
            A dictionary containing the OLS-regression values for
            * beta: Least-squares solution.
            * residuals: Sums of residuals; squared Euclidean 2-norm
            * rank: rank of the input matrix
            * singular_values: The singular values of input matrix
        monthly_avg: (pandas.Series)
            The climatological expectation for each calendar month of a year
        monthly_std: (pandas.Series)
            The standard deviation calculated from the observations and
            the climatological expectation for each calendard month of a year
        """
        # Resample the data to monthly means
        mu = ds[param].resample(time="M").mean()

        # Next, build the model
        ts = mu.values
        f = 1/12
        N = len(ts)
        t_in = np.arange(0, N, 1)
        t_out = t_in

        # Drop NaNs from the fit
        mask = np.isnan(ts)
        ts = ts[mask == False]
        t_in = t_in[mask == False]
        n = len(t_in)

        # Build the 2-cycle model
        X = [np.ones(n), np.sin(2*np.pi*f*t_in), np.cos(2*np.pi*f*t_in),
             np.sin(4*np.pi*f*t_in), np.cos(4*np.pi*f*t_in)]
        [beta, resid, rank, s] = np.linalg.lstsq(np.transpose(X), ts)
        self.regression = {
            "beta": beta,
            "residuals": resid,
            "rank": rank,
            "singular_values": s
        }

        # Calculate the two-cycle fitted data
        fitted_data = beta[0] + beta[1]*np.sin(2*np.pi*f*t_out) + beta[2]*np.cos(
            2*np.pi*f*t_out) + beta[3]*np.sin(4*np.pi*f*t_out) + beta[4]*np.cos(4*np.pi*f*t_out)
        fitted_data = pd.Series(fitted_data, index=mu.get_index("time"))
        self.fitted_data = fitted_data

        # Get the standard deviations
        self.calc_standard_deviations(ds, param)
