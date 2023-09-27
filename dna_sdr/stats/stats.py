"""
General statistical analysis functions

    Functions:
       __array_to_df(data, cond_lst, group_name, value_name) 
"""

import os
import itertools as it
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.stats.multicomp as mc
import statsmodels.stats.multitest as mt
from statsmodels.formula.api import ols


def __array_to_df(
    data: np.array,
    condition: list,
    group_name: str = "group",
    value_name: str = "value",
):
    """
    A hidden function that converts numpy array input into a pandas dataframe

    Parameters
    ----------
    data : np.array

    condition : list
        A list of condition presented
    group_name : str
        The name of the group array with default name of "group"
    value_name : str
        The name of the value array with default name of "value"

    Returns
    -------
    pd.DataFrame

    """

    if not condition:
        groups = list(range(1, len(data) + 1))
    else:
        groups = condition

    group_label = ""
    for index, group in enumerate(data):
        condition = [str(groups[index])]
        label = condition * len(group)
        if len(group_label) == 0:
            group_label = label
        else:
            group_label = group_label + label

    group_label_lst = list(group_label)
    comb_data = []
    for i in data:
        comb_data.extend(i)

    frame = pd.DataFrame({f"{group_name}": group_label_lst, f"{value_name}": comb_data})
    return frame, groups


def fence_rule(data: np.array):
    """
    xxx

    Parameters
    ----------
    data : np.array

    Returns
    -------
    fence_in : np.array

    fence_out_count : int

    """

    q_1, q_3 = np.quantile(data, [0.25, 0.75], method="median_unbiased")
    iqr = stats.iqr(data, interpolation="median_unbiased")
    upper = q_3 + 1.5 * iqr
    lower = q_1 - 1.5 * iqr
    fence_in = data[(lower <= data) & (upper >= data)]
    fence_out_count = len(data) - len(fence_in)
    return fence_in, fence_out_count


def modified_z_score(data: np.array):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    mad = stats.median_abs_deviation(data)
    med = np.median(data)
    mod_zs = abs(data - med) * 0.6745 / mad
    return mod_zs


def general_esd(data: np.array, poss_outlier_count: int = 5, alpha: float = 0.05):
    """
    xxx

    Parameters
    ----------
    data : np.array

    poss_outlier_count : int

    alpha : float

    Returns
    -------
    data : np.array

    """

    def test_stat(data):
        data_mean = np.mean(data)
        data_stdev = np.std(data, ddof=1)
        deviation = abs(data - data_mean)
        max_deviation = max(deviation)
        max_deviation_index = np.argmax(deviation)
        data_point = data[max_deviation_index]
        t_stats = max_deviation / data_stdev
        return t_stats, max_deviation_index, data_point

    def critical_value(data, alpha):
        sample_count = len(data)
        degree_freedom = sample_count - 2
        p_value = 1 - (alpha / (2 * sample_count))
        t_distribution = stats.t.ppf(p_value, degree_freedom)
        num = (sample_count - 1) * t_distribution
        deno = np.sqrt((degree_freedom + t_distribution**2) * sample_count)
        crit_value = num / deno
        return crit_value

    try:
        data = data.to_numpy()
    except AttributeError:
        pass

    index_lst = []
    max_index = 0
    data_copy = data.copy()
    data_point_lst = []
    test_stats_lst = []
    critical_value_lst = []

    for position in range(1, poss_outlier_count + 1):
        t_stat, position, point = test_stat(data_copy)
        crit = critical_value(data_copy, alpha)
        data_point_lst.append(point)
        test_stats_lst.append(t_stat)
        critical_value_lst.append(crit)

        if t_stat > crit:
            max_index = position

        data_copy = np.delete(data_copy, position)
        index_lst.append(position)

    result_table = pd.DataFrame(
        {
            "Data point": data_point_lst,
            "Test statistics (Ri)": test_stats_lst,
            "Critical value (λi)": critical_value_lst,
        }
    )
    print(result_table)
    print(
        f"There are {max_index} possible outliers in the data based on the general ESD test"
    )
    outlier_index = index_lst[0:max_index]

    for position in outlier_index:
        data = np.delete(data, position)

    return data


def outlier_detection(data: list, analysis_target: list):
    """
    xxx

    Parameters
    ----------
    data : list

    analysis_target: list

    Returns
    -------
    fence_lst : list

    post_gesd_lst : list

    """

    fence_lst = []
    post_gesd_lst = []
    for index, sample in enumerate(data):
        print(analysis_target[index])
        stats.probplot(sample, plot=plt, rvalue=True)
        plt.show()
        wstats, pvalue = stats.shapiro(sample)
        print(f"The W stats is {wstats:.4f}")
        p_value_printout(pvalue)
        print("")
        fence_in, out_count = fence_rule(sample)
        print(
            f"There are {out_count} possilbe outliers based on the fence rule (1.5*IQR)"
        )
        print("")
        if out_count != 0:
            wstats, pvalue = stats.shapiro(fence_in)
            print("Shapiro test after removal of potential outliers with fence rule")
            print(f"The W stats is {wstats:.4f}")
            p_value_printout(pvalue)
            print("")
        print("General ESD Test for outlier")
        data_post_gesd = general_esd(sample, 5, 0.05)
        print("")
        if not np.array_equiv(data_post_gesd, sample):
            wstats, pvalue = stats.shapiro(data_post_gesd)
            print("Shapiro test after removal of potential outliers with general ESD")
            print(f"The W stats is {wstats:.4f}")
            p_value_printout(pvalue)
            print("")
        fence_lst.append(fence_in)
        post_gesd_lst.append(data_post_gesd)

    return fence_lst, post_gesd_lst


def p_value_printout(p_value, alpha: float = 0.05):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    if p_value < alpha:
        print(
            f"The p value ({p_value:.3e}) is less than the defined level of significance ({alpha})"
        )
        print("The null hypothesis is rejected")
    else:
        print(
            f"The p value ({p_value:.3e}) is greater than the defined level of significance ({alpha})"
        )
        print("The null hypothesis is not rejected")
    return


def anova_test(
    data: pd.DataFrame,
    parameter: list,
    model: str,
    dependent: str,
    independent: str = "plate_loc",
):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    def anova_table(aov):
        aov["mean_sq"] = aov[:]["sum_sq"] / aov[:]["df"]

        aov["eta_sq"] = aov[:-1]["sum_sq"] / sum(aov["sum_sq"])

        aov["omega_sq"] = (
            aov[:-1]["sum_sq"] - (aov[:-1]["df"] * aov["mean_sq"][-1])
        ) / (sum(aov["sum_sq"]) + aov["mean_sq"][-1])

        cols = ["sum_sq", "df", "mean_sq", "F", "PR(>F)", "eta_sq", "omega_sq"]
        aov = aov[cols]
        return aov

    print("")
    model_fit = ols(model, data=data).fit()
    aov = sm.stats.anova_lm(model_fit, typ=2)
    aov_table = anova_table(aov)
    print(aov_table)

    print("")
    print("Shapiro test for normality for the ANOVA residual")
    wstats, pvalue = stats.shapiro(model_fit.resid)
    print(f"The W stats is {wstats:.4f}")
    p_value_printout(pvalue)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    stats.probplot(model_fit.resid, plot=ax, rvalue=True)
    ax.set_title("Probability plot of model's residuals", fontsize=20)
    plt.show()

    print("")
    print("Levene test for equal variance")
    levene = stats.levene(*parameter)
    print(f"Test statistics is {levene[0]:.4f}")
    p_value_printout(levene[-1])

    print("")
    if aov["PR(>F)"][0] < 0.05:
        comp = mc.MultiComparison(data[dependent], data[independent])
        post_hoc_res = comp.tukeyhsd()
        print(post_hoc_res.summary())
        print("")

    return


def kw_test(data, parameter: str, alpha: float = 0.05):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    print("Kruskal-Wallis H-test")
    _, p_value = stats.kruskal(*data)
    result = pd.DataFrame(
        {
            "test type": "Kruskal-Wallis H-test",
            "parameter": parameter,
            "p-value": f"{p_value:.3e}",
            "Reject Null?": p_value < alpha,
        },
        index=[0],
    )
    return result, p_value


def dunn_test(data, conditions: list, adj_method: str = "hs"):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    print("Dunn's Test")

    def omega(frame, group1, group2, tie_sum):
        total_obs = len(frame)
        data_1 = frame[frame["group"] == str(group1)]["value"]
        data_2 = frame[frame["group"] == str(group2)]["value"]
        ele_1 = (total_obs * (total_obs + 1)) / 12
        ele_2 = tie_sum / (12 * (total_obs - 1))
        ele_3 = 1 / len(data_1) + 1 / len(data_2)
        omega_value = np.sqrt((ele_1 - ele_2) * ele_3)
        return omega_value

    def mean_rank_diff(mean_rank_df, group1, group2):
        mean_rank1 = mean_rank_df[mean_rank_df["group"] == str(group1)]["rank"].iloc[0]
        mean_rank2 = mean_rank_df[mean_rank_df["group"] == str(group2)]["rank"].iloc[0]
        return mean_rank1 - mean_rank2

    frame, groups = __array_to_df(data, conditions)
    frame["rank"] = frame["value"].rank()
    frame_rank_mean = frame.groupby("group")["rank"].mean().reset_index()
    counts = frame.groupby("rank")["value"].agg("count").reset_index()
    tie_sum = np.sum(
        counts[counts["value"] != 1]["value"] ** 3
        - counts[counts["value"] != 1]["value"]
    )

    comparison_lst = list(it.combinations(groups, 2))
    p_value_lst = list()
    for comp_group in comparison_lst:
        omega_i = omega(frame, *comp_group, tie_sum)
        y_i = mean_rank_diff(frame_rank_mean, *comp_group)
        z_i = y_i / omega_i
        p_value = 2.0 * stats.norm.sf(np.abs(z_i))
        p_value_lst.append(p_value)

    hypo_test, p_adj, *_ = mt.multipletests(p_value_lst, method=adj_method)
    result_table = pd.DataFrame(
        {
            "Comparison group": comparison_lst,
            "Adjusted p-value": p_adj,
            "Reject Null?": hypo_test,
        }
    )

    result_table["Adjusted p-value"] = [
        f"{value:.3e}" for value in result_table["Adjusted p-value"]
    ]

    return result_table


def conover_iman_test(data, conditions: list, adj_method: str = "hs"):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    print("Conover-Iman test")

    def test_stats(frame, groups, h_corr, s2, group1, group2):
        rank_mean = frame.groupby("group")["rank"].mean().reset_index()
        mean_rank1 = rank_mean[rank_mean["group"] == str(group1)]["rank"].iloc[0]
        mean_rank2 = rank_mean[rank_mean["group"] == str(group2)]["rank"].iloc[0]
        mean_rank_diff = abs(mean_rank1 - mean_rank2)

        n1 = len(frame[frame["group"] == str(group1)])
        n2 = len(frame[frame["group"] == str(group2)])

        se = np.sqrt(
            s2
            * ((len(frame) - 1 - h_corr) / (len(frame) - len(groups)))
            * ((1 / n1) + (1 / n2))
        )
        tstats = mean_rank_diff / se
        return tstats

    h_corr, _ = stats.kruskal(*data)
    frame, groups = __array_to_df(data, conditions)
    total_obs = len(frame)
    frame["rank"] = frame["value"].rank()
    counts = frame.groupby("rank")["value"].agg("count").reset_index()

    tie_sum = np.sum(
        counts[counts["value"] != 1]["value"] ** 3
        - counts[counts["value"] != 1]["value"]
    )
    tie_correction = 1 - (tie_sum / (total_obs**3 - total_obs))

    if tie_correction == 1:
        s2 = (total_obs * (total_obs + 1)) / 12
    else:
        s2 = (1 / (total_obs - 1)) * (
            np.sum(frame["rank"] ** 2) - ((total_obs * ((total_obs + 1) ** 2)) / 4)
        )

    comparison_lst = list(it.combinations(groups, 2))
    p_value_lst = []
    for comp in comparison_lst:
        tstats = test_stats(frame, groups, h_corr, s2, *comp)
        p_value = 2.0 * stats.t.sf(np.abs(tstats), df=total_obs - len(groups))
        p_value_lst.append(p_value)

    hypo_test, p_adj, *_ = mt.multipletests(p_value_lst, method=adj_method)
    result_table = pd.DataFrame(
        {
            "Comparison group": comparison_lst,
            "Adjusted p-value": p_adj,
            "Reject Null?": hypo_test,
        }
    )

    result_table["Adjusted p-value"] = [
        f"{value:.3e}" for value in result_table["Adjusted p-value"]
    ]
    return result_table


def data_dist_plot(data: pd.DataFrame, mode: str):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    fig = plt.figure(figsize=(10, 10))

    if mode == "one_phase":
        ax1 = fig.add_subplot(121)
        sns.boxplot(
            x="plate_loc", y="rate", data=data, color="#99c2a2", width=0.2, ax=ax1
        )
        sns.swarmplot(x="plate_loc", y="rate", data=data, color="#7d0013", ax=ax1)

        ax2 = fig.add_subplot(122)
        sns.boxplot(
            x="plate_loc", y="plateau", data=data, color="#99c2a2", width=0.2, ax=ax2
        )
        sns.swarmplot(x="plate_loc", y="plateau", data=data, color="#7d0013", ax=ax2)

    else:
        sns.boxplot(x="plate_loc", y="k_rate", data=data, color="#99c2a2", width=0.2)
        sns.swarmplot(x="plate_loc", y="k_rate", data=data, color="#7d0013")

    plt.show()
    return


def parameter_array(data: pd.DataFrame, mode: str, analysis_target: list):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    if mode == "one_phase":
        rate_lst = list()
        plateau_lst = list()
        for target in analysis_target:
            target_df = data[data["plate_loc"] == target]
            rate_lst.append(target_df["rate"].to_numpy())
            plateau_lst.append(target_df["plateau"].to_numpy())
        return rate_lst, plateau_lst

    else:
        rate_lst = list()
        for target in analysis_target:
            target_df = data[data["plate_loc"] == target]
            rate_lst.append(target_df["k_rate"].to_numpy())
        return rate_lst


def main(analysis_target: list, mode: str):
    """
    xxx

    Parameters
    ----------
    data : np.array


    Returns
    -------
    mode_zs : float

    """

    fname = f"indiviual_Screen_{mode}_param.pkl"
    df = pd.read_pickle(fname)
    df.rename({"Condition": "plate_loc"}, axis=1, inplace=True)
    trig_info = pd.read_pickle("trig_info.pkl")
    trig_df = pd.merge(df, trig_info, on="plate_loc")

    trig_df.set_index("plate_loc", inplace=True)
    analysis_df = trig_df.loc[analysis_target]
    analysis_df.reset_index(inplace=True)
    return analysis_df


if __name__ == "__main__":
    # def anova_table(aov):
    #     aov["mean_sq"] = aov[:]["sum_sq"] / aov[:]["df"]

    #     aov["eta_sq"] = aov[:-1]["sum_sq"] / sum(aov["sum_sq"])

    #     aov["omega_sq"] = (
    #         aov[:-1]["sum_sq"] - (aov[:-1]["df"] * aov["mean_sq"][-1])
    #     ) / (sum(aov["sum_sq"]) + aov["mean_sq"][-1])

    #     cols = ["sum_sq", "df", "mean_sq", "F", "PR(>F)", "eta_sq", "omega_sq"]
    #     aov = aov[cols]
    #     return aov

    os.chdir("dna_sdr/pickles/")
    # df = pd.read_pickle("Screen_first_kinetic_param.pkl")
    # print(df)

    read_lst = ["P0_A1", "P2_A5", "P2_A6", "P0_A5", "P2_D6", "P2_D7"]
    MODE = "one_phase"
    same_loc = main(read_lst, MODE)
    loc_lst = []
    for loc in same_loc["mismatch_loc"]:
        loc_lst.append(loc[0])

    type_lst = []
    for types in same_loc["mismatch_type"]:
        type_lst.append(types[0])
    same_loc["mismatch_loc"] = loc_lst
    same_loc["mismatch_type"] = type_lst

    print(same_loc)

    # model = "plateau ~ C(mismatch_loc) + C(mismatch_type) + C(mismatch_loc):C(mismatch_type)"

    # model_fit = ols(model, data=same_loc).fit()
    # aov = sm.stats.anova_lm(model_fit, typ=2)
    # aov_table = anova_table(aov)
    # print(aov_table)
    # print("")
    # print("Shapiro test for normality for the ANOVA residual")
    # w, pvalue = stats.shapiro(model_fit.resid)
    # print(f"The W stats is {w:.4f}")
    # p_value_printout(pvalue)
    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(111)
    # normality_plot, stat = stats.probplot(model_fit.resid, plot=plt, rvalue=True)
    # ax.set_title("Probability plot of model's residuals", fontsize=20)
    # plt.show()

    # model = (
    #     "rate ~ C(mismatch_loc) + C(mismatch_type) + C(mismatch_loc):C(mismatch_type)"
    # )

    # model_fit = ols(model, data=same_loc).fit()
    # aov = sm.stats.anova_lm(model_fit, typ=2)
    # aov_table = anova_table(aov)
    # print(aov_table)
    # print("")
    # print("Shapiro test for normality for the ANOVA residual")
    # w, pvalue = stats.shapiro(model_fit.resid)
    # print(f"The W stats is {w:.4f}")
    # p_value_printout(pvalue)
    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(111)
    # normality_plot, stat = stats.probplot(model_fit.resid, plot=plt, rvalue=True)
    # ax.set_title("Probability plot of model's residuals", fontsize=20)
    # plt.show()

    # """
    # Graphical repersentation of the data distribution using boxplot and catagorical scatter plot
    # (swarmplot) using seaborn package.
    # """
    # data_dist_plot(same_loc)
    # parameter_lst = parameter_array(same_loc, mode, read_lst)
    # if len(parameter_lst) == 2:
    #     rate_lst, plateau_lst = parameter_lst

    # """
    # Outlier detection using fence rule (1.5*IQR) and generalized ESD test.
    # Shapiro test was also used to test the normality of the data set before and after
    # the removal of outliers. Generalized ESD test assumes the samples are approx.
    # normally distributed.
    # """
    # fence_in, post_gesd = outlier_detection(plateau_lst, read_lst)

    # """
    # ANOVA testing with normality plot, Shapiro test (normality test),
    # and Levene test (equal variance test) for assumpition testing.
    # """
    # model = "plateau ~ C(plate_loc)"
    # anova_test(same_loc, plateau_lst, model, "plateau")

    # """
    # Kruskal-Wallis H-test (non-parametric version of ANOVA) in case when the
    # ANOVA assumption was not satisified.
    # """
    # result, p_value = kw_test(plateau_lst, "plateau")
    # print(result)
    # if p_value < 0.05:
    #     print(dunn_test(plateau_lst, read_lst))
    #     print(conover_iman_test(plateau_lst, read_lst))
