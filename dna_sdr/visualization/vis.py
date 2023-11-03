"""
Module docstring
"""
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import pylab
import seaborn as sns
from dna_sdr.experimental.data_processing import time_list_generation
from dna_sdr.sequence import seq_utilits
from dna_sdr.sequence.dna_utilits import DNA

CB_color_cycle = [
    "#377eb8",
    "#ff7f00",
    "#4daf4a",
    "#f781bf",
    "#a65628",
    "#984ea3",
    "#999999",
    "#e41a1c",
    "#dede00",
]

params = {
    # "figure.figsize": [6.4, 4.8],
    "axes.labelsize": 24,
    "axes.titlesize": 16,
    "axes.linewidth": 2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelpad": 6,
    "axes.formatter.use_mathtext": True,
    "errorbar.capsize": 5,
    "xtick.labelsize": 20,
    "xtick.major.size": 6,
    "xtick.major.width": 2,
    "ytick.labelsize": 20,
    "ytick.major.size": 6,
    "ytick.major.width": 2,
}

pylab.rcParams.update(params)

marker_dict = {1: ".", 2: "x", 3: "s", 4: "|"}


def mismatched_maker_gen(comparision_dict: dict[int, int]) -> dict:
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    mismatched_loc = {k: v for k, v in comparision_dict.items() if v != 0}
    mismatched_maker_dict = {}

    for k, v in mismatched_loc.items():
        mismatched_maker_dict[k] = marker_dict[v]

    return mismatched_maker_dict


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


def seq_comp_plot(seq_1: str, seq_2: str, ax=None):
    """
    xxx

    Parameters
    ----------

    Returns
    -------

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


def swarm_box_plot(df: pd.DataFrame, x_target: str, param_type: str, ax=None):
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    if ax is None:
        ax = plt.gca()

    box_plot_info = {
        "x": x_target,
        "y": param_type,
        "data": df,
        "color": "#99c2a2",
        "width": 0.5,
        "linewidth": 1,
        "linecolor": "black",
        "flierprops": {"marker": "x"},
        "medianprops": {"linewidth": 1.3},
        "whiskerprops": {"linewidth": 1.5},
        "capprops": {"linewidth": 1.5},
    }
    swarm_plot_info = {
        "x": x_target,
        "y": param_type,
        "data": df,
        "color": "#7d0013",
        "size": 5,
    }

    ax = sns.boxplot(**box_plot_info)
    ax = sns.swarmplot(**swarm_plot_info)
    return ax


def release_grid_plot(
    data: pd.Series,
    result,
    max_value: float,
    ax=None,
    props: dict[str, str | float] | None = None,
):
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    if ax is None:
        ax = plt.gca()

    if props is None:
        props = {
            "boxstyle": "round",
            "facecolor": "wheat",
            "alpha": 0.5,
        }

    time = time_list_generation(60)
    time.insert(0, 0)
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
    PARAMS = f"./dna_sdr/pickles/individual_{TEST}_param.pkl"
    READ = ["P1_A1", "P1_A7", "P1_A8", "P1_A9", "P1_A10"]
    info_df: pd.DataFrame = pd.read_pickle(INFO)
    params_df: pd.DataFrame = pd.read_pickle(PARAMS)
    result_df = pd.merge(info_df, params_df, on="plate_loc")
    read_df = result_df[result_df["plate_loc"].isin(READ)]
    print(read_df)

    # fig = plt.figure(figsize=(10, 7.5))
    # ax1 = swarm_box_plot(read_df, "toehold", "rate")
    # ax1 = swarm_box_plot(read_df, "toehold", "k_rate")
    # ax1.set(xlabel="Number of toehold base (nt)", ylabel="Rate (nM/min$^{-1}$)")
    # plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    # plt.show()
    # fig.savefig(
    #     f"./dna_sdr/image/summery/{TEST}/toehold_kinetic_rate.pdf", format="pdf"
    # )
