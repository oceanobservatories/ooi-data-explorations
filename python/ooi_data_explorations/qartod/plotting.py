# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


# +
def plot_data_variable(ds, param, add_deployments=True):
    """Function to plot the timeseries with deployment info.
    
    Parameters
    ----------
    ds: (xarray.Dataset)
        An xarray dataset downloaded from OOINet
    param: (str)
        The parameter name of the data variable in the OOI
        dataset to plot
    add_deployments: (boolean)
        Also plot deployment information
        
    Returns
    -------
    fig, ax: (matplotlib figs)
        Figure and axis handles for the matplotlib image
    """
    
    # Initialize the plot
    fig, ax = plt.subplots(figsize=(12,8))
    
    # Plot the data
    if add_deployments:
        s = ax.scatter(ds["time"], ds[param], c=ds["deployment"])
    else:
        ax.scatter(ds["time"], ds[param], c="tab:blue")
        
    # Add in labels
    yavg, ystd = ds[param].mean(), ds[param].std()
    ax.grid()
    ax.set_ylabel(ds[param].attrs["long_name"])
    ax.set_title(ds.attrs["id"])
    
    # Add in the deployments with vertical lines with text for the deployments
    if add_deployments:
        deployments = np.unique(ds["deployment"])
        for depNum in deployments:
            dt = ds.where(ds["deployment"] == depNum, drop=True)["time"].min()
            ax.vlines(dt, yavg-5*ystd, yavg+5*ystd)
            ax.text(dt, yavg-4*ystd, str(depNum), fontsize=12)
        
        # Add in the legend
        m, l = s.legend_elements()[0], s.legend_elements()[1]
        l = ["Deployment " + x for x in l]
        legend = (m, l)
        ax.legend(*legend, edgecolor="black", loc="center left", bbox_to_anchor=(1, 0.5))
        
    return fig, ax
<<<<<<< HEAD
        
=======
>>>>>>> dosta_cgsn_updates

def plot_gross_range(ds, param, gross_range):
    """Plot the data with the associated climatology values.
    
    Parameters
    ----------
    ds: (xarray.Dataset)
        An xarray dataset downloaded from OOINet
    param: (str)
        The parameter name of the data variable in the OOI
        dataset to plot.
    gross_range: (qartod.gross_range object)
        An object containing the calculated gross_range values
        for the associated dataset and variable
    """
    # Initialize the data
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Get axis limits
    yavg, ystd = ds[param].mean(), ds[param].std()
    tmin, tmax = ds.time.min().values, ds.time.max().values

    # Scatter plot the data
    ax.plot(ds.time, ds[param], linestyle="", marker=".", color="tab:red")
    ax.fill_between([tmin, tmax], gross_range.suspect_min, gross_range.suspect_max, color="tab:red", alpha=0.3)
    ax.set_ylim(yavg-7*ystd, yavg+7*ystd)
    ax.grid()
    ax.set_ylabel(ds[param].attrs["long_name"])
    ax.set_title(ds.attrs["id"])
    fig.autofmt_xdate()
    
    return fig, ax
        
        
def plot_climatology(ds, param, climatology):
    """Plot the data with the associated climatology values.
    
    Parameters
    ----------
    ds: (xarray.Dataset)
        An xarray dataset downloaded from OOINet
    param: (str)
        The parameter name of the data variable in the OOI
        dataset to plot.
    climatology: (qartod.climatology object)
        An object containing the calculated climatology values
        for the associated dataset and variable
        """
    fig, ax = plt.subplots(figsize = (12, 8))
    
    # Observations
    ax.plot(ds.time, ds[param], marker=".", linestyle="", color="tab:red", zorder=0, label="Observations")
    yavg, ystd = np.mean(ds[param]), np.std(ds[param])
    ymin, ymax = yavg-ystd*5, yavg+ystd*5
    ax.set_ylim((ymin, ymax))
    
    # Standard Deviation +/- 3
    for t in climatology.fitted_data.index:
        t0 = pd.Timestamp(year=t.year, month=t.month, day=1)
        mu = climatology.monthly_fit.loc[t.month]
        std = climatology.monthly_std.loc[t.month]
        ax.hlines(mu, t0, t, color="black", linewidth=3, label="Climatological Fit")
        ax.fill_between([t0, t], [mu+3*std, mu+3*std], [mu-3*std, mu-3*std], color="tab:red", alpha=0.3, label="3*$\sigma$")
    ax.grid()
    
    # Add legend and labels
    handles, labels = ax.get_legend_handles_labels()[0][0:3], ax.get_legend_handles_labels()[1][0:3]
    ax.legend(handles, labels, fontsize=12)
    ax.set_title("-".join((ds.attrs["id"].split("-")[0:4])), fontsize=16)
    ax.set_ylabel(ds[param].attrs["long_name"], fontsize=16)
    fig.autofmt_xdate()
    
    return fig, ax
