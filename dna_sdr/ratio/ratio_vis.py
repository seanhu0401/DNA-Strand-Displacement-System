"""
This file is used to generate the graphs for the concentration studies using functions defined in `vis.py`
"""

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import pylab
from dna_sdr.visualization.vis import swarm_box_plot

params = {
    "figure.figsize": [9.18, 5],
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
    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Ratio_param.pkl"
    )

    T1 = "P0_A0"
    T2 = "P0_A2"

    read_df = params_df[
        (params_df["Trig1"].str.contains(T1)) & (params_df["Trig2"].str.contains(T2))
    ]
    read_df[["T1_Conc", "T2_Conc"]] = read_df["Concentration"].str.split(
        "_", expand=True
    )
    read_df[["T1_Conc", "T2_Conc"]] = read_df[["T1_Conc", "T2_Conc"]].apply(
        pd.to_numeric
    )

    one_phase_df = read_df[read_df["r_sq_curve"] >= 0.85]
    kinetic_df = read_df[read_df["r_sq_kin"] >= 0.85]

    pivot = pd.pivot_table(
        one_phase_df,
        values=["plateau", "rate"],
        index=["Concentration"],
        aggfunc=["mean", "std", "median", "count"],
    )

    # print(pivot)
    fig, axes = plt.subplots(2, 2, figsize=(9.18, 5), layout="constrained")
    swarm_box_plot(kinetic_df, "T1_Conc", "k_rate", ax=axes[0, 1])
    axes[0, 1].set(xlabel="Trigger 1 Ratio (%)", ylabel="Rate (nM/min)")
    swarm_box_plot(one_phase_df, "T1_Conc", "plateau", ax=axes[1, 0])
    swarm_box_plot(one_phase_df, "T1_Conc", "rate", ax=axes[1, 1])
    axes[1, 1].set(xlabel="Trigger 1 Ratio (%)", ylabel="Rate (min$^{-1}$)")
    axes[1, 0].set(xlabel="Trigger 1 Ratio (%)", ylabel="Quantity (nM)")
    plt.show()

    # fig.savefig(
    #     "./dna_sdr/image/summery/Ratio/T1_P0_A2_V3.pdf",
    #     format="pdf",
    # )
