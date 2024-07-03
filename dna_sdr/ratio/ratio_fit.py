import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
import pandas as pd
from pandas import DataFrame

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
pd.set_option("display.max_columns", 10)


def quantile_25(series: pd.Series):
    return series.quantile(0.25)


def quantile_75(series: pd.Series):
    return series.quantile(0.75)


def pred_plateau(
    parameter_df: DataFrame, conc_1: float, conc_2: float, trig_pair: tuple
):
    Trig1, Trig2 = trig_pair

    if Trig1 == "P0_A0":
        Trig1 = "T1"

    def exp_plateau(conc: float, plateau: float, k: float):
        return plateau * (1 - np.exp(-k * conc))

    trig1_pair = (parameter_df.loc[Trig1]["Plateau"], parameter_df.loc[Trig1]["Rate"])
    trig2_pair = (parameter_df.loc[Trig2]["Plateau"], parameter_df.loc[Trig2]["Rate"])

    ratio_plateau = exp_plateau(conc_1, *trig1_pair) + exp_plateau(conc_2, *trig2_pair)

    return ratio_plateau


def pred_rate(conc_1: float, rate_1: float, conc_2: float, rate_2: float):
    return


if __name__ == "__main__":
    CONC_RESULT = "./dna_sdr/pickles/conc_plateau_fit_result_unweighted_v4.pkl"
    RATE_RESULT = "./dna_sdr/pickles/conc_rate.pkl"

    fit_result_df = pd.read_pickle(CONC_RESULT)
    fit_result_df = fit_result_df.set_index("Trigger")
    print(fit_result_df)

    rate_df = pd.read_pickle(RATE_RESULT)
    mean_rate = rate_df["mean"]["rate"].rename("mean")
    std_rate = rate_df["std"]["rate"].rename("std")
    rates = pd.concat([mean_rate, std_rate], axis=1)

    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Ratio_param.pkl"
    )

    T1 = "P0_A0"
    T2 = "P0_A2"

    read_df = params_df[
        (params_df.Trig1.str.contains(T1)) & (params_df.Trig2.str.contains(T2))
    ]
    read_df[["T1_Conc", "T2_Conc"]] = read_df["Concentration"].str.split(
        "_", expand=True
    )
    read_df[["T1_Conc", "T2_Conc"]] = read_df[["T1_Conc", "T2_Conc"]].apply(
        pd.to_numeric
    )

    trig_pairs = (T1, T2)
    conc_pairs = list(set(zip(read_df.T1_Conc / 100, read_df.T2_Conc / 100)))
    conc_pairs.sort()

    pred_plateau_lst = []
    for conc_pair in conc_pairs:
        print(conc_pair)
        print(pred_plateau(fit_result_df, *conc_pair, trig_pair=trig_pairs))
        pred_plateau_lst.append(
            pred_plateau(fit_result_df, *conc_pair, trig_pair=trig_pairs)
        )

    one_phase_df = read_df[read_df["r_sq_curve"] >= 0.85]

    pivot = pd.pivot_table(
        one_phase_df,
        values=["rate", "plateau"],
        index=[
            "T1_Conc",
            "T2_Conc",
        ],
        aggfunc=["mean", "std", "median", quantile_25, quantile_75],
    )
    print(pivot)

    fig, axes = plt.subplots(
        2, 1, figsize=(4.6, 9.1), layout="constrained", sharex=True
    )
    types = list(set(one_phase_df["T1_Conc"]))
    types.sort()
    COUNT = 0

    axes[0].errorbar(
        types, pivot["mean"]["plateau"], pivot["std"]["plateau"], fmt="o", color="blue"
    )
    axes[1].errorbar(
        types, pivot["mean"]["rate"], pivot["std"]["rate"], fmt="o", color="blue"
    )

    axes[0].legend()
    axes[0].set(ylabel="Quantity (nM)")
    axes[1].set(xlabel="Trigger 1 Ratio (%)", ylabel="Rate (min$^{-1}$)")

    plt.show()
    # fig.savefig(
    #     f"./dna_sdr/image/summery/Ratio/.pdf",
    #     format="pdf",
    # )
