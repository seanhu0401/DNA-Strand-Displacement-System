"""
This file is used to generate the graphs for the concentration studies using functions defined in `vis.py`
"""

import pandas as pd
import warnings
import matplotlib.pyplot as plt
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
warnings.filterwarnings("ignore")

if __name__ == "__main__":
    READ = "T1"
    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Conc_param.pkl"
    )

    read_df = params_df[params_df["Trig"].str.contains(READ)]
    read_df["Concentration"] = read_df["Concentration"].astype(int)
    read_df["Concentration"] = read_df["Concentration"] / 100
    read_df = read_df[read_df["plateau"] < 700]
    curve_df = read_df[read_df["r_sq_curve"] >= 0.5]
    kinetic_df = read_df[read_df["r_sq_kin"] >= 0.5]

    fig, axes = plt.subplots(2, 2, figsize=(9.18, 5), layout="constrained")
    swarm_box_plot(kinetic_df, "Concentration", "k_rate", ax=axes[0, 1])
    axes[0, 1].set(xlabel="Trigger:Target Ratio (x)", ylabel="Rate (nM/min)")
    axes[0, 1].ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    swarm_box_plot(curve_df, "Concentration", "plateau", ax=axes[1, 0])
    swarm_box_plot(curve_df, "Concentration", "rate", ax=axes[1, 1])
    axes[1, 1].set(xlabel="Trigger:Target Ratio (x)", ylabel="Rate (min$^{-1}$)")
    axes[1, 0].set(xlabel="Trigger:Target Ratio (x)", ylabel="Quantity (nM)")
    plt.show()

    # fig.savefig(
    #     f"./dna_sdr/image/summery/Conc/{READ}_new_fit_method_V4.pdf",
    #     format="pdf",
    # )
