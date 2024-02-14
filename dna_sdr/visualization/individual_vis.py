"""
Module docstring
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import pylab

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

pylab.rcParams.update(params)


if __name__ == "__main__":
    time_lst = time_list_generation(60)
    MAX_VALUE = 500
    df: DataFrame = pd.read_pickle(
        "./dna_sdr/IO/Output/Individual/4WJ_HEX_Screen_P2_A5.pkl"
    )

    # y0 = [500, 500, 0, 0]
    # model = Model(sec_kin_fit, independent_vars=["t", "y0"])
    # params = model.make_params(k1=1e-07)

    fig, axes = plt.subplots(4, 3, figsize=(12, 8), layout="constrained")
    ROW = 0
    COLUMN = 0

    props = {
        "boxstyle": "round",
        "facecolor": "wheat",
        "alpha": 0.5,
    }

    for _, trial in df.items():
        conc = trial * 500
        if MAX_VALUE <= conc.max():
            MAX_VALUE = conc.max()
        y_axis_range = np.arange(0, MAX_VALUE + 50, 100)

        if conc.iloc[0] > 0:
            model = Model(one_phase_association)
            params = model.make_params(plateau=250, k=0.01, y_0=0)
        else:
            model = Model(lag_one_phase_association)
            params = model.make_params(plateau=250, k=0.01, time_0=100, y_0=0)

        result = model.fit(conc, params, time=time_lst)
        # result = model.fit(conc, params, t=time_lst, y0=y0)
        print(result.fit_report())

        axes[ROW, COLUMN].scatter(time_lst, conc)
        axes[ROW, COLUMN].plot(
            time_lst, result.best_fit, "-", label="best fit", color="#7d0013"
        )

        text_str = "\n".join(
            (
                f"$p.$: {result.best_values['plateau']:.2e}\xB1{result.params['plateau'].stderr:.2e}",
                f"$k$: {result.best_values['k']:.2e}\xB1{result.params['k'].stderr:.2e}",
                f"$r^2$ = {result.rsquared:.3f}",
            )
        )

        # text_str = "\n".join(
        #     (
        #         f"$k$: {result.best_values['k1']:.2e}\xB1{result.params['k1'].stderr:.2e}",
        #         f"$r^2$ = {result.rsquared:.3f}",
        #     )
        # )

        axes[ROW, COLUMN].text(
            0.50,
            0.45,
            text_str,
            transform=axes[ROW, COLUMN].transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=props,
        )

        if ROW % 3 == 0 and ROW != 0:
            COLUMN += 1
            ROW = 0
        else:
            ROW += 1

    fig.supylabel("Release quantity (nM)", fontsize=20)
    fig.supxlabel("Time (mins)", fontsize=20)
    # fig.suptitle(f"{name}", fontsize=24)
    plt.show()
