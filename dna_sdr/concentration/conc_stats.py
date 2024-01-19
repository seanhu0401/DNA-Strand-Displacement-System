"""
This file is used to perform stats analysis for the concentration studies using functions defined in `stats.py`
"""

import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import statsmodels.stats.multicomp as mc
from scipy import stats
from statsmodels.formula.api import ols

from dna_sdr.stats.stats import anova_table, p_value_printout

if __name__ == "__main__":
    READ = "A2"
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

    print("Plateau ANOVA")
    MODEL = "plateau ~ C(concentration)"
    model_fit = ols(MODEL, data=read_df).fit()
    aov = sm.stats.anova_lm(model_fit, typ=2)
    aov_table_plateau = anova_table(aov)
    print(aov_table_plateau)

    print("")
    print("Shapiro test for normality for the ANOVA residual")
    wstats, pvalue = stats.shapiro(model_fit.resid)
    print(f"The W stats is {wstats:.4f}")
    p_value_printout(pvalue)

    fig, axes = plt.subplots(2, 1, figsize=(9.18, 5), layout="constrained", sharex=True)
    normality_plot, stat = stats.probplot(model_fit.resid, plot=axes[0], rvalue=True)
    axes[0].set_title("Probability plot of model's residuals - Plateau", fontsize=20)

    comp = mc.MultiComparison(read_df["plateau"], read_df["concentration"])
    post_hoc_res = comp.tukeyhsd()
    tukey_post_hoc_plateau = pd.DataFrame(
        data=post_hoc_res._results_table.data[1:],
        columns=post_hoc_res._results_table.data[0],
    )
    print(tukey_post_hoc_plateau)

    print("Rate ANOVA")
    MODEL = "rate ~ C(concentration)"
    model_fit = ols(MODEL, data=read_df).fit()
    aov = sm.stats.anova_lm(model_fit, typ=2)
    aov_table_rate = anova_table(aov)
    print(aov_table_rate)

    print("")
    print("Shapiro test for normality for the ANOVA residual")
    wstats, pvalue = stats.shapiro(model_fit.resid)
    print(f"The W stats is {wstats:.4f}")
    p_value_printout(pvalue)

    normality_plot, stat = stats.probplot(model_fit.resid, plot=axes[1], rvalue=True)
    axes[1].set_title("Probability plot of model's residuals - Rate", fontsize=20)

    comp = mc.MultiComparison(read_df["rate"], read_df["concentration"])
    post_hoc_res = comp.tukeyhsd()
    tukey_post_hoc_rate = pd.DataFrame(
        data=post_hoc_res._results_table.data[1:],
        columns=post_hoc_res._results_table.data[0],
    )
    print(tukey_post_hoc_rate)

    plt.show()

    # MODEL = "k_rate ~ C(toehold)"
    # model_fit = ols(MODEL, data=df).fit()
    # aov = sm.stats.anova_lm(model_fit, typ=2)
    # aov_table_k_plateau = anova_table(aov)
    # print(aov_table_k_plateau)

    # print("")
    # print("Shapiro test for normality for the ANOVA residual")
    # wstats, pvalue = stats.shapiro(model_fit.resid)
    # print(f"The W stats is {wstats:.4f}")
    # p_value_printout(pvalue)

    # fig, axes = plt.subplots(2, 1, figsize=(9.18, 5), layout="constrained", sharex=True)
    # normality_plot, stat = stats.probplot(model_fit.resid, plot=axes[0], rvalue=True)
    # axes[0].set_title(
    #     "Probability plot of model's residuals - Rate constant",
    #     fontsize=20,
    # )

    # comp = mc.MultiComparison(df["k_rate"], df["toehold"])
    # post_hoc_res = comp.tukeyhsd()
    # tukey_post_hoc_k_rate = pd.DataFrame(
    #     data=post_hoc_res._results_table.data[1:],
    #     columns=post_hoc_res._results_table.data[0],
    # )
    # print(tukey_post_hoc_k_rate)
