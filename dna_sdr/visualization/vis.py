"""
Module docstring
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lmfit.model import ModelResult
from matplotlib import pylab
from matplotlib.axes import Axes

from dna_sdr.experimental.data_processing import time_list_generation
from dna_sdr.sequence import seq_utilits
from dna_sdr.sequence.dna_utilits import DNA

params = {
    # "figure.figsize": [6.4, 4.8],
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "axes.linewidth": 1,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelpad": 4,
    "axes.formatter.use_mathtext": True,
    "errorbar.capsize": 5,
    "xtick.labelsize": 12,
    "xtick.major.size": 4,
    "xtick.major.width": 1,
    "ytick.labelsize": 12,
    "ytick.major.size": 4,
    "ytick.major.width": 1,
}

pylab.rcParams.update(params)

marker_dict = {1: ".", 2: "x", 3: "s", 4: "|"}


def mismatched_maker_gen(comparison_dict: dict[int, int]) -> dict[int, str]:
    """The function `mismatched_maker_gen` takes a dictionary of integers and returns a new dictionary
    where the keys are the same as the original dictionary, but the values are obtained by looking up
    the corresponding value in another dictionary called `marker_dict`.

    Parameters
    ----------
    comparison_dict : dict[int, int]
        A dictionary where the keys are integers representing locations and the values are integers
    representing markers.

    Returns
    -------
        a dictionary called `mismatched_maker_dict`.

    """
    mismatched_loc = {k: v for k, v in comparison_dict.items() if v != 0}
    mismatched_maker_dict = {}

    for k, v in mismatched_loc.items():
        mismatched_maker_dict[k] = marker_dict[v]

    return mismatched_maker_dict


def pt_to_inch(x_pt: float, y_pt: float) -> tuple[float, float]:
    """The function `pt_to_inch` converts a given point value to inches.

    Parameters
    ----------
    x_pt : float
        The parameter `x_pt` represents the value in points for the x-coordinate. Points are a unit of
    measurement commonly used in typography and graphic design, where 1 point is equal to 1/72 of an
    inch.
    y_pt : float
        The parameter `y_pt` represents the y-coordinate in points. Points are a unit of measurement
    commonly used in typography and graphic design.

    Returns
    -------
        a tuple of two floats, representing the values of x_inches and y_inches.

    """
    x_inches = x_pt / 72
    y_inches = y_pt / 72
    return x_inches, y_inches


def seq_comp_plot(seq_1: str, seq_2: str, ax: Axes | None = None) -> Axes:
    """The `seq_comp_plot` function takes two DNA sequences as input and plots a comparison between them,
    highlighting any mismatches.

    Parameters
    ----------
    seq_1 : str
        The `seq_1` parameter is a string representing the first sequence for comparison.
    seq_2 : str
        The `seq_2` parameter is a string representing the second sequence that you want to compare with
    `seq_1`.
    ax : Axes | None
        The `ax` parameter is an optional parameter of type `Axes` or `None`. It represents the matplotlib
    Axes object on which the sequence comparison plot will be drawn. If `ax` is `None`, a new Axes
    object will be created using `plt.gca()`.

    Returns
    -------
        the Axes object that is used to create the plot.

    """

    if ax is None:
        ax = plt.gca()

    comp_dict = seq_utilits.sequence_comparison(seq_1, seq_2)
    seq_tuple = (DNA(seq_1), DNA(seq_2))
    base_loc_lst = list(range(seq_tuple[0].base_count))
    y_temp = [1] * len(base_loc_lst)
    ax.plot(base_loc_lst, y_temp)

    mismatched_maker = mismatched_maker_gen(comp_dict)

    marker_lst = [v for _, v in marker_dict.items()]
    for marker in marker_lst:
        mark_dict = {k: v for k, v in mismatched_maker.items() if v == marker}
        if bool(mark_dict):
            loc_lst = [loc for loc in mark_dict if mark_dict[loc] == marker]
            ax.plot(base_loc_lst, y_temp, marker, markevery=loc_lst, markersize=12)

    ax.set_xlabel("Base number")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax


def scatter_plot(
    release_df: pd.DataFrame,
    test_type: str,
    label: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """The scatter_plot function creates a scatter plot with error bars using data from a DataFrame.

    Parameters
    ----------
    release_df : pd.DataFrame
        A pandas DataFrame containing the data for the scatter plot. It should have at least two columns,
    with the second column being the x-axis values for the scatter plot.
    test_type : str
        The `test_type` parameter is a string that specifies the type of test data to be plotted. It is
    used to filter the columns of the `release_df` DataFrame to select the relevant data for the scatter
    plot.
    label : str | None
        The `label` parameter is used to provide a label for the scatter plot. It is an optional parameter
    and can be set to `None` if no label is needed.
    ax : Axes | None
        The `ax` parameter is an optional parameter of type `Axes`. It represents the matplotlib Axes
    object on which the scatter plot will be drawn. If no `ax` object is provided, the function will
    create a new Axes object using `plt.gca()` (get current axes) and use that for

    Returns
    -------
        an instance of the `Axes` class from the `matplotlib.pyplot` module.

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


def regression_plot(
    release_df: pd.DataFrame, model, fit_params, label, ax: Axes | None = None
) -> Axes:
    """The function `regression_plot` plots the fitted regression line on a given axis using the provided
    model and fit parameters.

    Parameters
    ----------
    release_df
        The release_df parameter is a pandas DataFrame that contains the data for the regression plot. It
    should have at least two columns, where the first column represents the x-axis values and the second
    column represents the y-axis values.
    model
        The model is a function that represents the regression model you want to plot. It takes the x-axis
    values as the first argument and any additional parameters required for the model. The function
    should return the predicted y-values based on the input x-values and model parameters.
    fit_params
        The `fit_params` parameter is a list of parameters that are used to fit the regression model. These
    parameters are specific to the chosen regression model and are used to estimate the relationship
    between the independent variable(s) and the dependent variable.
    label
        The label parameter is a string that represents the label for the plotted line. It is used to
    identify the line in the legend of the plot.
    ax
        The ax parameter is an optional parameter that represents the matplotlib axes object on which the
    regression plot will be plotted. If no ax object is provided, the function will create a new axes
    object using plt.gca().

    Returns
    -------
        the matplotlib axes object `ax`.

    """
    if ax is not None:
        pass
    else:
        ax = plt.gca()

    x_axis = release_df[release_df.columns.values[1]]
    fitted = model(x_axis, *fit_params)

    ax.plot(x_axis, fitted, label=label)

    return ax


def swarm_box_plot(
    df: pd.DataFrame,
    x_target: str,
    param_type: str,
    ax: Axes | None = None,
    hue: str | None = None,
    palette: tuple | None = None,
    legend: bool = True,
) -> Axes:
    """The `swarm_box_plot` function creates a combined box plot and swarm plot using the Seaborn library
    in Python.

    Parameters
    ----------
    df : pd.DataFrame
        The `df` parameter is a pandas DataFrame that contains the data for the box plot and swarm plot.
    x_target : str
        The x_target parameter is the name of the column in the DataFrame that you want to use as the
    x-axis variable for the box plot and swarm plot.
    param_type : str
        The parameter `param_type` represents the variable or parameter that you want to visualize in the
    box plot and swarm plot. It could be any numerical variable in your DataFrame `df` that you want to
    compare across different categories or groups defined by the `x_target` variable.
    ax : Axes | None
        The `ax` parameter is an optional parameter that represents the matplotlib Axes object on which the
    box plot and swarm plot will be drawn. If no `ax` object is provided, the function will create a new
    Axes object using `plt.gca()` (get current axes) and use that for plotting.

    Returns
    -------
        the Axes object that contains the box plot and swarm plot.

    """

    if ax is None:
        ax = plt.gca()

    box_plot_info = {
        "x": x_target,
        "y": param_type,
        "data": df,
        "color": "#99c2a2",
        "width": 0.5,
        "linewidth": 0.5,
        "linecolor": "black",
        "flierprops": {"marker": "x"},
        "medianprops": {"linewidth": 0.8},
        "whiskerprops": {"linewidth": 1.0},
        "capprops": {"linewidth": 0.8},
        "ax": ax,
        "fill": True,
        "hue": hue,
        "legend": legend,
    }
    swarm_plot_info = {
        "x": x_target,
        "y": param_type,
        "data": df,
        "color": "#7d0013",
        "size": 3,
        "ax": ax,
        "legend": False,
        "dodge": True,
        "hue": hue,
    }
    if palette is None:
        sns.boxplot(**box_plot_info)
        sns.swarmplot(**swarm_plot_info)
    else:
        sns.boxplot(palette=sns.color_palette(palette[0]), **box_plot_info)
        sns.swarmplot(palette=sns.color_palette(palette[1]), **swarm_plot_info)
    return ax


def release_grid_plot(
    data: pd.Series,
    result: ModelResult,
    max_value: float,
    time: list[float] | None = None,
    ax: Axes | None = None,
    props: dict[str, str | float] | None = None,
) -> Axes:
    """The `release_grid_plot` function plots data points and a model evaluation on a grid, with additional
    information displayed in a text box.

    Parameters
    ----------
    data : pd.Series
        The `data` parameter is a pandas Series object that contains the data points to be plotted.
    result : ModelResult
        The `result` parameter is an instance of the `ModelResult` class. It contains the result of a model
    fitting operation, such as the best-fit values of the model parameters (`result.best_values`), the
    covariance matrix (`result.covar`), and the coefficient of determination (`result.rs
    max_value : float
        The `max_value` parameter represents the maximum value for the y-axis range in the grid plot. It is
    used to set the upper limit of the y-axis range in the plot.
    time : list[float] | None
        The `time` parameter is a list of float values representing the time points for the data.
    ax : Axes | None
        The `ax` parameter is an optional argument of type `Axes` from the `matplotlib.axes` module. It
    represents the axes on which the grid plot will be drawn. If no `ax` object is provided, the
    function will use the current axes (`plt.gca()`) as the default value
    props : dict[str, str | float] | None
        The `props` parameter is a dictionary that contains properties for the text box that is displayed
    on the plot. It can have the following keys:

    Returns
    -------
        The function `release_grid_plot` returns an instance of the `Axes` class from the
    `matplotlib.pyplot` module.

    """
    if ax is None:
        ax = plt.gca()

    if props is None:
        props = {
            "boxstyle": "round",
            "facecolor": "wheat",
            "alpha": 0.5,
        }

    if time is None:
        time = time_list_generation(60)

    y_axis_range = np.arange(0, max_value + 50, 100)

    text_str = "\n".join(
        (
            f"$k_1$: {result.best_values['k1']:.2e}\xB1{np.sqrt(result.covar[0][0]):.2e}",
            f"$r^2$ = {result.rsquared:.3f}",
            f"Max={data.max():.1f}",
        )
    )
    ax.plot(time, result.eval(), "--r")
    ax.scatter(time, data, marker=".")
    ax.set_ylim(0, max_value + 50)
    ax.set_yticks(y_axis_range)
    ax.text(
        0.50,
        0.45,
        text_str,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=props,
    )

    return ax


if __name__ == "__main__":
    TEST = "Screen"
    INFO = "./dna_sdr/pickles/trig_info.pkl"
    PARAMS = f"./dna_sdr/pickles/individual_{TEST}_params.pkl"
    ### Mismatch type 6/7
    # READ = ["P2_A3", "P2_A7", "P2_D4", "P2_D5"]
    ###
    # READ = ["P0_A1", "P2_A2", "P2_A4", "P2_D3", "P2_D7", "P2_D11"]
    ###
    # READ = [
    #     "P0_A1",
    #     "P0_A5",
    #     "P2_A2",
    #     "P2_A4",
    #     "P2_A5",
    #     "P2_A6",
    #     "P2_A10",
    #     "P2_D3",
    #     "P2_D6",
    #     "P2_D7",
    #     "P2_D11",
    # ]

    # sns.scatterplot(
    #     x="plateau",
    #     y="rate",
    #     data=curve_df,
    #     hue="mismatch_type",
    #     ax=axes[0, 0],
    #     palette="bright",
    #     style="mismatch_loc",
    # )
    # axes[0, 0].set_yscale("log")

    ### Mismatch location variation plot
    # fig, axes = plt.subplots(2, 2, figsize=(9.18, 5), layout="constrained")
    # swarm_box_plot(curve_df, "mismatch_loc", "plateau", ax=axes[1, 0])
    # swarm_box_plot(curve_df, "mismatch_loc", "rate", ax=axes[1, 1])
    # swarm_box_plot(kinetic_df, "mismatch_loc", "k_rate", ax=axes[0, 1])
    # axes[1, 1].set(xlabel=x_label, ylabel="Rate (min$^{-1}$)")
    # axes[1, 0].set(xlabel=x_label, ylabel="Quantity (nM)")
    # axes[1, 1].set_yscale("log")
    # axes[0, 1].set(xlabel=x_label, ylabel="Rate (nM/min)")
    # axes[0, 1].set_yscale("log")
    # fig.savefig(f"./dna_sdr/image/summery/{TEST}/mismatch_loc_V1.pdf", format="pdf")
    # fig.savefig(
    #     f"./dna_sdr/image/summery/{TEST}/mismatch_loc_mismatch_6_7.pdf", format="pdf"
    # )

    plt.show()
