"""
This file is used to generate the graphs for the screening studies using functions defined in `vis.py`
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import pylab
import matplotlib.ticker as mticker
import seaborn as sns

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
pd.set_option("display.max_columns", 20)
pd.options.display.float_format = "{:.3e}".format
formatter = mticker.ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-3, 2))


if __name__ == "__main__":
    # Toehold variation
    # READ = ["P1_A1", "P1_A7", "P1_A8", "P1_A9", "P1_A10"]
    # Mismatch type 1-3 same location (17)
    READ = ["P0_A1", "P2_A5", "P2_A6"]
    # Mismatch type 4-6 same location (18)
    # READ = ["P2_A7", "P3_A3", "P3_A4"]
    # Mismatch type 7-9 same location (15)
    # READ = ["P2_A3", "P3_A5", "P3_A6"]
    # Mismatch type 1 varying locations
    # READ = ["P0_A1", "P2_A2", "P2_A4", "P2_D3", "P2_D7", "P2_D11"]
    # Mismatch type 2/8 varying location
    # READ = ["P0_A5", "P2_A1", "P2_A5", "P2_A10", "P2_D9", "P3_A1", "P3_A5"]
    # Mismatch group 1 at two location (17, 29)
    # READ = ["P0_A1", "P2_A5", "P2_A6", "P0_A5", "P2_D6", "P2_D7"]
    # Mismatch group 1 at location 29
    # READ = ["P0_A5", "P2_D6", "P2_D7"]

    INFO = "./dna_sdr/pickles/trig_info.pkl"
    info_df: pd.DataFrame = pd.read_pickle(INFO)
    info_df = info_df.rename({"plate_loc": "Trig"}, axis=1)
    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Screen_param.pkl"
    )
    result_df = pd.merge(info_df, params_df, on="Trig")
    read_df = result_df[result_df["Trig"].isin(READ)]

    # toehold = "toehold"
    # toehold_label = "Toehold (nt)"

    mis_type = "mismatch_type"
    mis_type_label = "Mismatch Type"

    mis_loc = "mismatch_loc"
    mis_loc_label = "Mismatch Location"

    type_lst = [int(i[0]) for i in read_df["mismatch_type"]]
    read_df["mismatch_type"] = type_lst
    loc_lst = [int(i[0]) for i in read_df["mismatch_loc"]]
    read_df["mismatch_loc"] = loc_lst

    curve_df = read_df[read_df["r_sq_curve"] >= 0.5]
    kinetic_df = read_df[read_df["r_sq_kin"] >= 0.5]

    # Toehold length variation plot
    # fig, axes = plt.subplots(1, 3, figsize=(9.18, 5), layout="constrained")
    # swarm_box_plot(kinetic_df, toehold, "k_rate", ax=axes[0])
    # swarm_box_plot(curve_df, toehold, "plateau", ax=axes[1])
    # swarm_box_plot(curve_df, toehold, "rate", ax=axes[2])
    # axes[0].set(xlabel=toehold_label, ylabel="Rate (nM/min)")
    # axes[0].set_yscale("log")
    # axes[1].set(xlabel=toehold_label, ylabel="Quantity (nM)")
    # axes[2].set(xlabel=toehold_label, ylabel="Rate (min$^{-1}$)")
    # axes[2].set_yscale("log")
    # fig.savefig("./dna_sdr/image/summery/Screen/toehold_V3.pdf", format="pdf")

    fig, axes = plt.subplots(2, 2, figsize=(9.18, 5), layout="constrained")

    # Mismatch type same location plot
    swarm_box_plot(kinetic_df, mis_type, "k_rate", ax=axes[0, 1])
    swarm_box_plot(curve_df, mis_type, "plateau", ax=axes[1, 0])
    swarm_box_plot(curve_df, mis_type, "rate", ax=axes[1, 1])
    axes[0, 1].set(xlabel=mis_type_label, ylabel="Rate (nM/min)")
    axes[0, 1].yaxis.set_major_formatter(formatter)
    axes[0, 1].set_yscale("log")
    axes[1, 1].set(xlabel=mis_type_label, ylabel="Rate (min$^{-1}$)")
    axes[1, 1].set_yscale("log")
    axes[1, 0].set(xlabel=mis_type_label, ylabel="Quantity (nM)")
    # fig.savefig(
    #     "./dna_sdr/image/summery/screen/mismatch_type_loc17_V2.pdf", format="pdf"
    # )

    # Mismatch location variation plot
    # swarm_box_plot(kinetic_df, mis_type, "k_rate", ax=axes[0, 1])
    # swarm_box_plot(curve_df, mis_type, "plateau", ax=axes[1, 0])
    # swarm_box_plot(curve_df, mis_type, "rate", ax=axes[1, 1])
    # axes[0, 1].set(xlabel=mis_type_label, ylabel="Rate (nM/min)")
    # axes[0, 1].set_yscale("log")
    # axes[1, 1].set(xlabel=mis_type_label, ylabel="Rate (min$^{-1}$)")
    # axes[1, 1].set_yscale("log")
    # axes[1, 0].set(xlabel=mis_type_label, ylabel="Quantity (nM)")
    # fig.savefig(
    #     "./dna_sdr/image/summery/screen/mismatch_loc_type2_8_V3.pdf", format="pdf"
    # )

    # 1-3 Mismatch Type - Two different location
    # palette = ("pastel", "dark")
    # swarm_box_plot(
    #     kinetic_df, mis_type, "k_rate", ax=axes[0, 1], hue=mis_loc, palette=palette, legend=False
    # )
    # swarm_box_plot(
    #     curve_df, mis_type, "plateau", ax=axes[1, 0], hue=mis_loc, palette=palette
    # )
    # swarm_box_plot(
    #     curve_df, mis_type, "rate", ax=axes[1, 1], hue=mis_loc, palette=palette, legend=False
    # )
    # axes[0, 1].set(xlabel=mis_type_label, ylabel="Rate (nM/min)")
    # axes[0, 1].set_yscale("log")
    # axes[1, 1].set(xlabel=mis_type_label, ylabel="Rate (min$^{-1}$)")
    # axes[1, 1].set_yscale("log")
    # axes[1, 0].set(xlabel=mis_type_label, ylabel="Quantity (nM)")
    # sns.move_legend(
    #     axes[1, 0],
    #     "lower center",
    #     bbox_to_anchor=(0.5, 1),
    #     ncol=4,
    #     title=None,
    #     frameon=False,
    # )
    # fig.savefig(
    #     "./dna_sdr/image/summery/screen/mismatch_type_1_3_V2.pdf", format="pdf"
    # )

    plt.show()
