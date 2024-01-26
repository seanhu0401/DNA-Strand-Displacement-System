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
    RESULT = "./dna_sdr/pickles/conc_plateau_fit_result.pkl"
    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Conc_param.pkl"
    )

    params_df["Concentration"] = params_df["Concentration"].astype(int)
    params_df["Concentration"] = params_df["Concentration"] / 100
    params_df = params_df[params_df["plateau"] < 700]

    pivot_rate = pd.pivot_table(
        params_df,
        values=["rate", "k_rate"],
        index=["Trig"],
        aggfunc=["mean", "std", "count"],
    )
    pivot_rate.to_pickle("./dna_sdr/pickles/conc_rate.pkl")

    pivot_plateau = pd.pivot_table(
        params_df,
        values=["plateau"],
        index=[
            "Trig",
            "Concentration",
        ],
        aggfunc=["mean", "std", "var", "count"],
    )

    fig, axes = plt.subplots(3, 1, figsize=(4.6, 9.1), layout="constrained", sharex=True)
    types = ["T1", "P0_A1", "P0_A2"]
    COUNT = 0

    # text_file = open(FNAME, "a")
    result_list = []

    for trig_type in types:
        mean = pivot_plateau.loc[trig_type]["mean"].squeeze()
        std = pivot_plateau.loc[trig_type]["std"].squeeze()
        weight = 1 / pivot_plateau.loc[trig_type]["var"].squeeze()
        concentration = list(pivot_plateau.loc[trig_type].index)
        params_list = [trig_type]

        model = Model(exp_plateau)
        params = model.make_params(plateau=250, k=0.01)

        result = model.fit(mean, params, conc=concentration)

        # text_file.write(trig_type + "\n")
        for name, par in result.params.items():
            params_list.extend([par.value, par.stderr])
        #     text_file.write(f"{name} = {par.value:.3f} +/- {par.stderr:.3f}\n")
        # text_file.write(f"r_squared = {result.rsquared:.5f}\n")
        result_list.append(params_list)

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
    conc_plateau_fit_df = pd.DataFrame(result_list, columns=["Trigger", "Plateau", "SE_Plateau", "Rate", "SE_Rate"])
    conc_plateau_fit_df.to_pickle(RESULT)
