"""
Module docstring
"""

from typing import Sequence
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import pylab
from matplotlib.axes import Axes

from lmfit import Model
from pandas import DataFrame

from dna_sdr.experimental.data_processing import time_list_generation
from dna_sdr.fitting.fit import (
    one_phase_association,
    lag_one_phase_association,
    sec_kin_fit,
)


params = {
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


props = {
    "boxstyle": "round",
    "facecolor": "wheat",
    "alpha": 0.5,
}


pylab.rcParams.update(params)


def individual_fitting_figure(
    y: pd.Series,
    fit_type: str,
    ax: Axes | None = None,
    y0: Sequence[float | int] | None = None,
):
    if ax is None:
        ax = plt.gca()

    if fit_type == "one_phase":
        if y.iloc[0] > 0:
            model = Model(one_phase_association)
            params = model.make_params(plateau=250, k=0.01, y_0=0)
        else:
            model = Model(lag_one_phase_association)
            params = model.make_params(plateau=250, k=0.01, time_0=100, y_0=0)
        result = model.fit(y, params, time=time_lst)
        text_str = "\n".join(
            (
                f"$p.$: {result.best_values['plateau']:.2f}\xB1{result.params['plateau'].stderr:.2f}",
                f"$k$: {result.best_values['k']:.2e}\xB1{result.params['k'].stderr:.2e}",
                f"$y_0$: {result.best_values['y_0']:.2f}\xB1{result.params['y_0'].stderr:.2f}",
                f"$r^2$ = {result.rsquared:.3f}",
            )
        )

    elif fit_type == "second_kinetic":
        model = Model(sec_kin_fit, independent_vars=["t", "y0"])
        params = model.make_params(k1=1e-07)
        result = model.fit(y, params, t=time_lst, y0=y0)
        text_str = "\n".join(
            (
                f"$k$: {result.best_values['k1']:.2e}\xB1{result.params['k1'].stderr:.2e}",
                f"$r^2$ = {result.rsquared:.3f}",
            )
        )

    else:
        raise ValueError("Unknown fit type")

    ax.scatter(time_lst, y)
    ax.plot(time_lst, result.best_fit, "-", label="best fit", color="#7d0013")

    ax.text(
        0.50,
        0.55,
        text_str,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=props,
    )

    return


if __name__ == "__main__":
    CONC = "100"
    TARGET = "T1"
    time_lst = time_list_generation(60)
    MAX_VALUE = 500
    df: DataFrame = pd.read_pickle(
        f"./dna_sdr/IO/Output/Individual/4WJ_HEX_Conc_{TARGET}_{CONC}.pkl"
    )
    chunks = math.ceil(len(df.columns) / 12)

    y0 = [500, 500, 0, 0]

    ROW = 0
    COLUMN = 0

    FIT = "one_phase"
    n = len(df.columns) // chunks

    list_df = [df.iloc[:, i : i + n] for i in range(0, len(df.columns), n)]

    for figure in range(chunks):
        fig, axes = plt.subplots(4, 3, figsize=(12, 8), layout="constrained")
        for _, trial in list_df[figure].items():
            conc = trial * 500
            if MAX_VALUE <= conc.max():
                MAX_VALUE = conc.max()
            y_axis_range = np.arange(0, MAX_VALUE + 50, 100)

            individual_fitting_figure(conc, FIT, axes[ROW, COLUMN], y0)

            if ROW % 3 == 0 and ROW != 0:
                COLUMN += 1
                ROW = 0
            else:
                ROW += 1

        COLUMN = 0
        ROW = 0

        fig.supylabel("Release quantity (nM)", fontsize=20)
        fig.supxlabel("Time (mins)", fontsize=20)

        plt.show()

        fig.savefig(
            f"./dna_sdr/image/individual/{TARGET}_{CONC}_Conc_{FIT}_figure_p{figure+1}.pdf",
            format="pdf",
        )
