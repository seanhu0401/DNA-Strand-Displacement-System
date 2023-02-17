import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {
    "axes.labelsize": 20,
    "axes.spines.top": False,
    "axes.spines.right": True,
    "errorbar.capsize": 5,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
}
pylab.rcParams.update(params)


def pt_to_inch(x_pt, y_pt):
    x_inches = x_pt / 72
    y_inches = y_pt / 72
    return x_inches, y_inches


def scatter_plot(release_df, type, label=None, ax=None):
    if ax is None:
        ax = plt.gca()

    x_axis = release_df[release_df.columns.values[1]]
    plateau = release_df.filter(regex=(type))
    mean = plateau.filter(regex=("mean")).squeeze()
    std = plateau.filter(regex=("std")).squeeze()

    ax.scatter(x_axis, mean, label=label, s=50)
    ax.errorbar(x_axis, mean, std, ls="none", elinewidth=2)

    return ax


def regression_plot(release_df, model, params, label, ax=None):
    if ax is None:
        ax = plt.gca()

    x_axis = release_df[release_df.columns.values[1]]
    fitted = model(x_axis, *params)

    ax.plot(x_axis, fitted, label=label)

    return ax
