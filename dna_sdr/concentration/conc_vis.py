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
    READ = "T1"
    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Conc_param.pkl"
    )

    read_df = params_df[params_df["plate_loc"].str.contains(READ)]
    read_df[["type", "concentration"]] = read_df["plate_loc"].str.split(
        "_", expand=True
    )
    read_df["concentration"] = read_df["concentration"].astype(int)
    read_df["concentration"] = read_df["concentration"] / 100
    read_df = read_df[read_df["plateau"] < 700]

    fig, axes = plt.subplots(2, 2, figsize=(9.18, 5), layout="constrained")
    swarm_box_plot(read_df, "concentration", "k_rate", ax=axes[0, 1])
    axes[0, 1].set(xlabel="Trigger:Target Ratio (x)", ylabel="Rate (nM/min)")
    swarm_box_plot(read_df, "concentration", "plateau", ax=axes[1, 0])
    swarm_box_plot(read_df, "concentration", "rate", ax=axes[1, 1])
    axes[1, 1].set(xlabel="Trigger:Target Ratio (x)", ylabel="Rate (min$^{-1}$)")
    axes[1, 0].set(xlabel="Trigger:Target Ratio (x)", ylabel="Quantity (nM)")
    plt.show()

    # fig.savefig(
    #     f"./dna_sdr/image/summery/Conc/P0_{READ}_V1.pdf",
    #     format="pdf",
    # )
