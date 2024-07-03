"""
This file is used to perform fitting for the results of the concentration study
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
import pandas as pd
from lmfit import Model

from dna_sdr.fitting.model_function import one_phase_association


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

props = {
    "boxstyle": "round",
    "facecolor": "wheat",
    "alpha": 0.5,
}


if __name__ == "__main__":
    FNAME = "./dna_sdr/result/Conc/conc_plateau_fit_result_unweighted_v4.csv"
    RESULT = "./dna_sdr/pickles/conc_plateau_fit_result_unweighted_v4.pkl"
    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Conc_param.pkl"
    )

    types = ["T1", "P0_A1", "P0_A2"]

    params_df["Concentration"] = params_df["Concentration"].astype(int)
    params_df["Concentration"] = params_df["Concentration"] / 100
    params_df = params_df[
        (params_df["plateau"] < 700) & (params_df["r_sq_curve"] > 0.5)
    ]

    # params_df['Trig'] = pd.Categorical(params_df.Trig, ordered=True, categories=types)
    # params_df = params_df.sort_values("Trig")

    pivot_rate = pd.pivot_table(
        params_df,
        values=["rate", "k_rate"],
        index=["Trig", "Concentration"],
        aggfunc=["mean", "std", "count"],
    )
    # pivot_rate.to_pickle("./dna_sdr/pickles/conc_rate.pkl")
    # pivot_rate.to_csv("./dna_sdr/result/Conc/conc_rate_v1.csv")
    # print(pivot_rate.reindex(index=types))

    pivot_plateau = pd.pivot_table(
        params_df,
        values=["plateau"],
        index=["Trig", "Concentration"],
        aggfunc=["mean", "std", "var", "count"],
    )

    # fig, axes = plt.subplots(
    #     3, 1, figsize=(4.6, 9.1), layout="constrained", sharex=True
    # )
    # COUNT = 0

    # result_list = []

    # for trig_type in types:
    #     mean = pivot_plateau.loc[trig_type]["mean"].squeeze()
    #     std = pivot_plateau.loc[trig_type]["std"].squeeze()
    #     weight = 1 / pivot_plateau.loc[trig_type]["var"].squeeze()
    #     concentration = list(pivot_plateau.loc[trig_type].index)
    #     params_list = [trig_type]

    #     model = Model(one_phase_association)
    #     # params = model.make_params(plateau=250, k=0.01, y_0=dict(value=70, min=0))
    #     params = model.make_params(plateau=250, k=5)
    #     # params = model.make_params(plateau=250, k=0.01)

    #     # result = model.fit(mean, params, conc=concentration, weights=1 / std**2)
    #     result = model.fit(mean, params, conc=concentration)
    #     print(result.fit_report())

    #     for name, par in result.params.items():
    #         params_list.extend([par.value, par.stderr])
    #     result_list.append(params_list)

    #     text_str = "\n".join(
    #         (
    #             f"$p.$: {result.best_values['plateau']:.2f}\xB1{result.params['plateau'].stderr:.2f}",
    #             f"$k$: {result.best_values['k']:.2f}\xB1{result.params['k'].stderr:.3f}",
    #             f"$r^2$ = {result.rsquared:.3f}",
    #         )
    #     )

    #     axes[COUNT].errorbar(
    #         concentration, mean, std, fmt="o", color="blue", label="experimental"
    #     )
    #     axes[COUNT].plot(
    #         concentration, result.best_fit, "-", label="best fit", color="#7d0013"
    #     )

    #     axes[COUNT].text(
    #         0.50,
    #         0.45,
    #         text_str,
    #         transform=axes[COUNT].transAxes,
    #         fontsize=10,
    #         verticalalignment="top",
    #         bbox=props,
    #     )

    #     COUNT += 1

    # axes[0].legend()
    # axes[0].set(ylabel="Quantity (nM)")
    # axes[1].set(ylabel="Quantity (nM)")
    # axes[2].set(xlabel="Trigger:Target Ratio (x)", ylabel="Quantity (nM)")
    # plt.show()
    # fig.savefig(
    #     "./dna_sdr/image/summery/Conc/plateau_fit_unweighted_V4.pdf",
    #     format="pdf",
    # )
    # conc_plateau_fit_df = pd.DataFrame(
    #     result_list,
    #     # columns=["Trigger", "Plateau", "SE_Plateau", "Rate", "SE_Rate", "Yo", "SE_Yo"],
    #     columns=["Trigger", "Plateau", "SE_Plateau", "Rate", "SE_Rate"],
    # )
    # conc_plateau_fit_df = pd.DataFrame(
    #     result_list, columns=["Trigger", "Plateau", "SE_Plateau", "Rate", "SE_Rate"]
    # )
    # print(conc_plateau_fit_df)
    # conc_plateau_fit_df.to_pickle(RESULT)
    # conc_plateau_fit_df.to_csv(FNAME)

    fig, axes = plt.subplots(
        3, 1, figsize=(4.6, 9.1), layout="constrained", sharex=True
    )
    COUNT = 0

    result_list = []

    for trig_type in types:
        mean = pivot_rate.loc[trig_type]["mean"].squeeze()["rate"]
        std = pivot_rate.loc[trig_type]["std"].squeeze()["rate"]
        concentration = list(pivot_rate.loc[trig_type].index)
        params_list = [trig_type]

        axes[COUNT].errorbar(
            concentration, mean, std, fmt="o", color="blue", label="experimental"
        )

        COUNT += 1

    axes[0].legend()
    axes[0].set(ylabel="Rate (min$^{-1}$)")
    axes[1].set(ylabel="Rate (min$^{-1}$)")
    axes[2].set(xlabel="Trigger:Target Ratio (x)", ylabel="Rate (min$^{-1}$)")
    plt.show()
    fig.savefig("./dna_sdr/image/summery/Conc/rate_V1.pdf", format="pdf")
