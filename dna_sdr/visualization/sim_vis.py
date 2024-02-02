"""
Module Docstring
"""

import matplotlib.pyplot as plt
from matplotlib import pylab

params = {
    "axes.labelsize": 20,
    "axes.spines.top": False,
    "axes.spines.right": True,
    "errorbar.capsize": 5,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
}
pylab.rcParams.update(params)


def pt_to_inch(x_pt: float, y_pt: float):
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    x_inches = x_pt / 72
    y_inches = y_pt / 72
    return x_inches, y_inches


def scatter_plot(release_df, test_type, label=None, ax=None):
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    if ax is None:
        ax = plt.gca()

    x_axis = release_df[release_df.columns.values[1]]
    plateau = release_df.filter(regex=test_type)
    mean = plateau.filter(regex="mean").squeeze()
    std = plateau.filter(regex="std").squeeze()

    ax.scatter(x_axis, mean, label=label, s=50)
    ax.errorbar(x_axis, mean, std, ls="none", elinewidth=2)

    return ax


def regression_plot(release_df, model, fit_params, label, ax=None):
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    if ax is None:
        ax = plt.gca()

    x_axis = release_df[release_df.columns.values[1]]
    fitted = model(x_axis, *fit_params)

    ax.plot(x_axis, fitted, label=label)

    return ax


# TODO: Add graphing function that incorporates scatter and regression plot into one plot
