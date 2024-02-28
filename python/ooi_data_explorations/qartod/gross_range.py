import os
import numpy as np
import pandas as pd
import xarray as xr
import dask
from dask.diagnostics import ProgressBar
from scipy.stats import normaltest


class GrossRange():
    """Gross Range fitting process for QARTOD.

    For a given parameter in a dataset, calculate the gross range QARTOD values
    for the data, and fromat the data into a qcConfig object optimized for use
    with Axiom-implemented QARTOD gross range test.

    Example
    -------
    from qartod.gross_range import GrossRange
    gross_range = GrossRange(fail_min=200, fail_max=1200)
    gross_range.fit(pco2_dataset, "pco2_seawater")
    gross_range.make_qcConfig()
    gross_range.qcConfig = {'qartod':
        {'gross_range_test':
            {'suspect_span': [200, 767.5], 'fail_span': [200, 1200]}}}
    """

    def __init__(self, fail_min, fail_max):
        """Init the Gross Range with the relevant fail min/max.

        Parameters
        ----------
        fail_min: (float)
            The minimum value for the given dataset parameter below which the
            data is fail.
        fail_max: (float)
            The maximum value for the given dataset parameter above which the
            data is fail.
        """
        self.fail_min = fail_min
        self.fail_max = fail_max

    def fit(self, ds, param, sigma=3, check_normality=False, **kwargs):
        """Fit suspect range with specified standard deviation.

        Parameters
        ----------
        ds: (xarray.DataSet)
            An xarray datasets containing the given data variable to be fitted,
            with a primary dimension of time.
        param: (string)
            The name of the data variable from the given dataset to fit
        sigma: (float: 3)
            The number of standard deviations for calculating the suspect range
        check_normality: (boolean: False)
            Check if the data is distributed normally. If not, utilize percentiles
            instead of mean and standard deviation.
        """
        # First, get the applicable data from the dataset
        da = ds[param]
        
        # Next, filter out data which falls outside of the fail ranges
        da = self.filter_fail_range(da)
        
        normal_flag = True
        # Next, check to see if data is normally distributed
        if check_normality:
            print(f"----- Testing {param} data for normality -----")
            pnorm = self.check_normality(da)
            if pnorm < 0.01:
                normal_flag = False
                print(f"----- {param} is not normally distributed -----")
            else:
                pass
        
        # If data is normally distributed, use mean and standard deviation
        if normal_flag:
            # Calculate the mean and standard deviation
            avg = np.nanmean(da)
            std = np.nanstd(da)

            # Calculate the suspect range
            suspect_min = avg-sigma*std
            suspect_max = avg+sigma*std
            source = f"User range based on the mean +- {sigma} standard deviations of all observations."
        else:
            # Even with a log-normal transformation, the data is not normally distributed, so we will
            # set the user range using percentiles that approximate the Empirical Rule
            if sigma == 3:
                lower_p = 0.15
                upper_p = 100-lower_p
            elif sigma == 4:
                lower_p = 3.167E-3
                upper_p = 100-lower_p
            elif sigma == 5:
                lower_p = 2.867E-5
                upper_p = 100-lower_p
            else:
                pass
            print(f"Using percentiles from {lower_p} to {upper_p}")
            # Calculate the suspect min/max based on percentiles equivalent to standard deviations
            suspect_min = np.nanpercentile(da, lower_p)
            suspect_max = np.nanpercentile(da, upper_p)
            p_range = np.round(100 - lower_p*2, 1)
            source = f"User range based on percentiles of the observations, which are not normally distributed. Percentiles were chosen to cover {p_range}% of the data, approximating the Empirical Rule."
            

        # If the suspect ranges are outside the fail ranges, set
        # suspect ranges to the fail_ranges
        if suspect_min < self.fail_min:
            suspect_min = self.fail_min
        if suspect_max > self.fail_max:
            suspect_max = self.fail_max

        # Save the results
        self.suspect_min = np.round(suspect_min, decimals=5)
        self.suspect_max = np.round(suspect_max, decimals=5)
        self.source = source

    def filter_fail_range(self, da):
        """Filter out values which fall outside the fail range."""
        mask = np.where((da > self.fail_min) & (da < self.fail_max) & (~np.isnan(da)))[0]
        da = da[mask]
        return da

    def check_normality(self, da, n_choices=5000, n_iters=1000):
        """Check if the data is normally distributed"""
        random_choice = dask.delayed(np.random.choice)
        vals = []
        for i in range(n_iters):
            vals.append(random_choice(da, n_choices))
        # Now compute
        with ProgressBar():
            pvals = dask.compute(*vals)
        # Take the mean of the pvalues of the normaltest
        pnorm = [normaltest(v).pvalue for v in pvals]
        return np.mean(pnorm)

    def make_qcConfig(self):
        """Build properly formatted qcConfig object for qartod gross range."""
        self.qcConfig = {
            "qartod": {
                "gross_range_test": {
                    "suspect_span": [self.suspect_min, self.suspect_max],
                    "fail_span": [self.fail_min, self.fail_max]
                }
            }
        }
