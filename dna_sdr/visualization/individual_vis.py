"""
Module docstring
"""

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import pylab
from matplotlib.axes import Axes

from pandas import DataFrame

from dna_sdr.experimental.data_processing import time_list_generation
from dna_sdr.fitting.fit import one_phase_fit, kinetic_fit


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
time_lst = time_list_generation(60)


def individual_fitting_figure(
    y: pd.Series,
    fit_type: str,
    test: str | None = None,
    condition: str | None = None,
    ax: Axes | None = None,
) -> None:
    if ax is None:
        ax = plt.gca()

    if fit_type == "one_phase":
        result = one_phase_fit(y)
        text_str = "\n".join(
            (
                f"$p.$: {result.best_values['plateau']:.2f}",
                f"$k$: {result.best_values['k']:.2e}",
                f"$y_0$: {result.best_values['y_0']:.2f}",
                f"$r^2$ = {result.rsquared:.3f}",
            )
        )

    elif fit_type == "second_kinetic":
        result = kinetic_fit(y, test, condition)
        # print(result.fit_report())
        try:
            text_str = "\n".join(
                (
                    f"$k$: {result.best_values['k1']:.2e}\xB1{result.params['k1'].stderr:.2e}",
                    f"$r^2$ = {result.rsquared:.3f}",
                )
            )
        except TypeError:
            text_str = "\n".join(
                (
                    f"$k$: {result.best_values['k1']:.2e}",
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
    conc_lst = [
        "",
        # "25",
        # "50",
        # "75",
        # "100",
        # "125",
        # "150",
        # "175",
        # "200",
        # "300",
        # "400",
        # "500",
        # "0_100",
        # "25_75",
        # "50_50",
        # "75_25",
        # "100_0",
    ]

    # TARGET = "P0_A0_P0_A2"
    TARGET = "P2_B9"
    TEST = "Screen"

    MAX_VALUE = 500
    FIT = "one_phase"
    # FIT = "second_kinetic"

    for CONC in conc_lst:
        print(CONC)
        df: DataFrame = pd.read_pickle(
            f"./dna_sdr/IO/Output/Individual/4WJ_HEX_{TEST}_{TARGET}.pkl"
            # f"./dna_sdr/IO/Output/Individual/4WJ_TR_{TEST}_{TARGET}.pkl"
            # f"./dna_sdr/IO/Output/Individual/4WJ_HEX_{TEST}_{TARGET}_{CONC}.pkl"
        )

        chunks = math.ceil(len(df.columns) / 12)

        ROW = 0
        COLUMN = 0

        n = len(df.columns) // chunks

        list_df = [df.iloc[:, i : i + n] for i in range(0, len(df.columns), n)]

        for figure in range(chunks):
            fig, axes = plt.subplots(4, 3, figsize=(12, 8), layout="constrained")
            for _, trial in list_df[figure].items():
                conc = trial * 500
                if MAX_VALUE <= conc.max():
                    MAX_VALUE = conc.max()
                y_axis_range = np.arange(0, MAX_VALUE + 50, 100)

                individual_fitting_figure(conc, FIT, TEST, CONC, axes[ROW, COLUMN])

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

            # fig.savefig("sample_release.pdf", format="pdf", transparent=True)

            fig.savefig(
                f"./dna_sdr/image/individual/{TARGET}_Screen_{FIT}_figure_p{figure+1}.pdf",
                format="pdf",
            )

            # fig.savefig(
            #     f"./dna_sdr/image/individual/{TARGET}_{CONC}_Conc_{FIT}_figure_p{figure+1}.pdf",
            #     format="pdf",
            # )
