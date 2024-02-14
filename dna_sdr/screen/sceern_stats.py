import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib import pylab
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.weightstats import ttest_ind
from pingouin import welch_anova, pairwise_gameshowell


pd.set_option("display.max_columns", 10)

params = {
    "axes.spines.top": False,
    "axes.spines.right": False,
}

pylab.rcParams.update(params)


def analysis_dataframe_prep(
    target: str | list[str], input_df: pd.DataFrame
) -> pd.DataFrame:
    """Prepare the dataframe for further data and stats analysis.

    Parameters
    ----------
    target: str | list[str]
        The `target` parameter represents trigger name or a list of trigger name interested.
    input_df: pd.DataFrame
        The `input_df` parameter represents the initial parameter dataframe object.

    Returns
    -------
        The function `analysis_dataframe_prep` returns an instance of `DataFrame` class from
        the `pandas` module.

    """
    INFO = "./dna_sdr/pickles/trig_info.pkl"
    info_df: pd.DataFrame = pd.read_pickle(INFO)
    info_df = info_df.rename({"plate_loc": "Trig"}, axis=1)
    result_df = pd.merge(info_df, input_df, on="Trig")
    target_df = result_df[result_df["Trig"].isin(target)]
    return target_df


def stats_analysis_routine(
    model: str,
    data: pd.DataFrame,
    independent: str,
    dependent: str,
    ax1: Axes | None = None,
    ax2: Axes | None = None,
) -> tuple:

    if ax1 is None:
        ax1 = plt.gca()
    if ax2 is None:
        ax2 = plt.gca()

    model_fit = ols(model, data=data).fit()

    welch = welch_anova(data, dv=dependent, between=independent)
    print(welch)

    games_howell = pairwise_gameshowell(data, dv=dependent, between=independent)
    print(games_howell)

    print(stats.shapiro(model_fit.resid))

    # Show the ANOVA residual against normal distribution
    stats.probplot(model_fit.resid, plot=ax1, rvalue=True)
    ax1.set_title(f"Probability plot of model's residuals - {dependent}", fontsize=14)

    ax2.scatter(data[independent], model_fit.resid)
    ax2.set_title(f"Variance plot of model's residuals - {dependent}", fontsize=14)

    return welch, games_howell


def cohens_d(sample1, sample2) -> float:
    return (np.mean(sample1) - np.mean(sample2)) / np.sqrt(
        (np.var(sample1) + np.var(sample2)) / 2  # type: ignore
    )


if __name__ == "__main__":
    # Toehold variation
    # READ = ["P1_A1", "P1_A7", "P1_A8", "P1_A9", "P1_A10"]
    # Mismatch type 1-3 same location (17)
    # READ = ["P0_A1", "P2_A5", "P2_A6"]
    # Mismatch type 4-6 same location (18)
    # READ = ["P2_A7", "P3_A3", "P3_A4"]
    # Mismatch type 7-9 same location (15)
    # READ = ["P2_A3", "P3_A5", "P3_A6"]
    # Mismatch type 1-3 same location (29)
    # READ = ["P0_A5", "P2_D6", "P2_D7"]
    # Mismatch type 1 varying locations
    # READ = ["P0_A1", "P2_A2", "P2_A4", "P2_D3", "P2_D7", "P2_D11"]
    # Mismatch type 2/8 varying location
    READ = ["P0_A5", "P2_A1", "P2_A5", "P2_A10", "P2_D9", "P3_A1", "P3_A5"]

    INFO = "./dna_sdr/pickles/trig_info.pkl"
    info_df: pd.DataFrame = pd.read_pickle(INFO)
    info_df = info_df.rename({"plate_loc": "Trig"}, axis=1)
    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Screen_param.pkl"
    )
    read_df = analysis_dataframe_prep(READ, params_df)
    read_df["rate"] = np.log(read_df["rate"])
    read_df["k_rate"] = np.log(read_df["k_rate"])

    mis_type = "mismatch_type"
    mis_type_label = "Mismatch Type"
    mis_loc = "mismatch_loc"
    mis_loc_label = "Mismatch Location"

    type_lst = [i[0] for i in read_df["mismatch_type"]]
    read_df["mismatch_type"] = type_lst
    loc_lst = [i[0] for i in read_df["mismatch_loc"]]
    read_df["mismatch_loc"] = loc_lst

    curve_df = read_df[read_df["r_sq_curve"] >= 0.5]
    kinetic_df = read_df[read_df["r_sq_kin"] >= 0.5]

    fig, axes = plt.subplots(3, 1, figsize=(6, 9.18), layout="constrained", sharex=True)
    fig2, axes2 = plt.subplots(
        3, 1, figsize=(6, 9.18), layout="constrained", sharex=True
    )

    # One way ANOVA toehold length variation
    # models = ["k_rate ~ C(toehold)", "plateau ~ C(toehold)", "rate ~ C(toehold)"]
    # aov_table_k_rate, games_howell_k_rate = stats_analysis_routine(
    #     models[0], kinetic_df, "toehold", "k_rate", axes[0], axes2[0]
    # )
    # aov_table_plateau, games_howell_plateau = stats_analysis_routine(
    #     models[1], curve_df, "toehold", "plateau", axes[1], axes2[1]
    # )
    # aov_table_rate, games_howell_rate = stats_analysis_routine(
    #     models[2], curve_df, "toehold", "rate", axes[2], axes2[2]
    # )

    # Mismatch type same location
    # models = [
    #     "k_rate ~ C(mismatch_type)",
    #     "plateau ~ C(mismatch_type)",
    #     "rate ~ C(mismatch_type)",
    # ]
    # aov_table_k_rate, games_howell_k_rate = stats_analysis_routine(
    #     models[0], kinetic_df, "mismatch_type", "k_rate", axes[0], axes2[0]
    # )
    # aov_table_plateau, games_howell_plateau = stats_analysis_routine(
    #     models[1], curve_df, "mismatch_type", "plateau", axes[1], axes2[1]
    # )
    # aov_table_rate, games_howell_rate = stats_analysis_routine(
    #     models[2], curve_df, "mismatch_type", "rate", axes[2], axes2[2]
    # )

    # Same mismatch type varying locations
    # models = [
    #     "k_rate ~ C(mismatch_loc)",
    #     "plateau ~ C(mismatch_loc)",
    #     "rate ~ C(mismatch_loc)",
    # ]
    # aov_table_k_rate, games_howell_k_rate = stats_analysis_routine(
    #     models[0], kinetic_df, "mismatch_loc", "k_rate", axes[0], axes2[0]
    # )
    # aov_table_plateau, games_howell_plateau = stats_analysis_routine(
    #     models[1], curve_df, "mismatch_loc", "plateau", axes[1], axes2[1]
    # )
    # aov_table_rate, games_howell_rate = stats_analysis_routine(
    #     models[2], curve_df, "mismatch_loc", "rate", axes[2], axes2[2]
    # )

    # 2 sample t-test w/ unequal variance
    # type_2 = curve_df[curve_df["mismatch_type"] == 2]
    # type_8 = curve_df[curve_df["mismatch_type"] == 8]
    # k_rate_ttest = list(ttest_ind(type_2["k_rate"], type_8["k_rate"], usevar="unequal"))
    # plateau_ttest = list(
    #     ttest_ind(type_2["plateau"], type_8["plateau"], usevar="unequal")
    # )
    # rate_ttest = list(ttest_ind(type_2["rate"], type_8["rate"], usevar="unequal"))
    # ttest_df = pd.DataFrame(
    #     [k_rate_ttest, plateau_ttest, rate_ttest],
    #     columns=["tstats", "pvalue", "df"],
    # )
    # print(cohens_d(type_2["k_rate"], type_8["k_rate"]))
    # print(cohens_d(type_2["plateau"], type_8["plateau"]))
    # print(cohens_d(type_2["rate"], type_8["rate"]))
    # ttest_df.to_csv("./dna_sdr/result/type_2_8_t_test.csv")

    # plt.show()

    # EXPORT = "mismatch_2_8_vary_location"
    # fig.savefig(
    #     f"./dna_sdr/image/summery/Screen/{EXPORT}_prob_plot_v3.pdf",
    #     format="pdf",
    # )
    # fig2.savefig(
    #     f"./dna_sdr/image/summery/Screen/{EXPORT}_variance_plot_v1.pdf",
    #     format="pdf",
    # )
    # aov_table = pd.concat([aov_table_k_rate, aov_table_plateau, aov_table_rate])
    # aov_table.to_csv(f"./dna_sdr/result/{EXPORT}_v3.csv")
    # mc_table = pd.concat([games_howell_k_rate, games_howell_plateau, games_howell_rate])
    # mc_table.to_csv(f"./dna_sdr/result/{EXPORT}_gameshowell_v1.csv")
