"""
This file is used to perform fitting for the results of the concentration study
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
import pandas as pd
from lmfit import Model

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


def exp_plateau(conc: float, plateau: float, k: float):
    return plateau * (1 - np.exp(-k * conc))


if __name__ == "__main__":
    FNAME = "./dna_sdr/result/Conc/conc_plateau_fit_result.txt"
    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Conc_param.pkl"
    )

    params_df[["type", "concentration"]] = params_df["plate_loc"].str.split(
        "_", expand=True
    )
    params_df["concentration"] = params_df["concentration"].astype(int)
    params_df["concentration"] = params_df["concentration"] / 100
    params_df = params_df[params_df["plateau"] < 700]

    pivot_rate = pd.pivot_table(
        params_df,
        values=["rate", "k_rate"],
        index=["type"],
        aggfunc=["mean", "std", "count"],
    )

    pivot_plateau = pd.pivot_table(
        params_df,
        values=["plateau"],
        index=[
            "type",
            "concentration",
        ],
        aggfunc=["mean", "std", "var", "count"],
    )

    fig, axes = plt.subplots(3, 1, figsize=(4.6, 9.1), layout="constrained", sharex=True)
    types = ["T1", "A1", "A2"]
    COUNT = 0
    # text_file = open(FNAME, "a")

    for trig_type in types:
        mean = pivot_plateau.loc[trig_type]["mean"].squeeze()
        std = pivot_plateau.loc[trig_type]["std"].squeeze()
        weight = 1 / pivot_plateau.loc[trig_type]["var"].squeeze()
        concentration = list(pivot_plateau.loc[trig_type].index)

        model = Model(exp_plateau)
        params = model.make_params(plateau=250, k=0.01)

        result = model.fit(mean, params, conc=concentration)

        # text_file.write(trig_type + "\n")
        # for name, par in result.params.items():
        #     text_file.write(f"{name} = {par.value:.3f} +/- {par.stderr:.3f}\n")
        # text_file.write(f"r_squared = {result.rsquared:.5f}\n")

        axes[COUNT].errorbar(concentration, mean, std, fmt="o", color="blue", label="experimental")
        axes[COUNT].plot(concentration, result.best_fit, "-", label="best fit", color="#7d0013")

        COUNT += 1

    # text_file.close()
    axes[0].legend()
    axes[0].set(ylabel="Quantity (nM)")
    axes[1].set(ylabel="Quantity (nM)")
    axes[2].set(xlabel="Trigger:Target Ratio (x)", ylabel="Quantity (nM)")
    # plt.show()
    # fig.savefig(
    #     f"./dna_sdr/image/summery/Conc/plateau_fit_unweighted_V1.pdf",
    #     format="pdf",
    # )
