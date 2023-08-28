import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.stats.multicomp as mc
import statsmodels.stats.multitest as mt
from statsmodels.formula.api import ols
import scipy.stats as stats
import itertools as it


def __array_to_df(data: np.array, group_name: str = "group", value_name: str = "value"):
    groups = list(range(1, len(data) + 1))
    group_label = ""
    for group in range(len(data)):
        label = str(groups[group]) * len(data[group])
        if len(group_label) == 0:
            group_label = label
        else:
            group_label = group_label + label

    group_label_lst = list(group_label)
    comb_data = list()
    for i in data:
        comb_data.extend(i)

    frame = pd.DataFrame(
        {"{}".format(group_name): group_label_lst, "{}".format(value_name): comb_data}
    )
    return frame, groups


def fence_rule(data: np.array):
    Q1, Q3 = np.quantile(data, [0.25, 0.75], method="median_unbiased")
    IQR = stats.iqr(data, interpolation="median_unbiased")
    upper = Q3 + 1.5 * IQR
    lower = Q1 - 1.5 * IQR
    fence_in = data[(lower <= data) & (upper >= data)]
    fence_out_count = len(data) - len(fence_in)
    return fence_in, fence_out_count


def general_ESD(data, poss_outlier_count: int = 2, alpha: float = 0.05):
    def test_stat(data):
        data_mean = np.mean(data)
        data_stdev = np.std(data, ddof=1)
        deviation = abs(data - data_mean)
        max_deviation = max(deviation)
        max_deviation_index = np.argmax(deviation)
        print(
            "The current data point tested is {:.3f}".format(data[max_deviation_index])
        )
        stats = max_deviation / data_stdev
        return stats, max_deviation_index

    def critical_value(data, alpha, poss_outlier_count):
        n = len(data)
        df = n - 2
        p = 1 - (alpha / (2 * n))
        t_distribution = stats.t.ppf(p, df)
        num = (n - 1) * t_distribution
        deno = np.sqrt((df + t_distribution**2) * n)
        crit_value = num / deno
        return crit_value

    try:
        data = data.to_numpy()
    except AttributeError:
        pass

    index_lst = list()
    max_index = 0
    data_copy = data.copy()

    for iter in range(1, poss_outlier_count + 1):
        print("Test {}".format(iter))
        stat, index = test_stat(data_copy)
        crit = critical_value(data_copy, alpha, poss_outlier_count)
        print("Test statistics (R{}): {:.3f}".format(iter, stat))
        print("Critical value (λ{}): {:.3f}".format(iter, crit))
        print()

        if stat > crit:
            max_index = iter

        data_copy = np.delete(data_copy, index)
        index_lst.append(index)

    print("There are {} possible outliers in the data".format(max_index))
    outlier_index = index_lst[0:max_index]

    for index in outlier_index:
        data = np.delete(data, index)

    return data


def p_value_printout(p_value, alpha: float = 0.05):
    if p_value < alpha:
        print(
            "The p value ({:.3e}) is less than the defined level of significance ({})".format(
                p_value, alpha
            )
        )
        print("The null hypothesis is rejected")
    else:
        print(
            "The p value ({:.3e}) is greater than the defined level of significance ({})".format(
                p_value, alpha
            )
        )
        print("The null hypothesis is not rejected")
    return


def anova_table(aov):
    aov["mean_sq"] = aov[:]["sum_sq"] / aov[:]["df"]

    aov["eta_sq"] = aov[:-1]["sum_sq"] / sum(aov["sum_sq"])

    aov["omega_sq"] = (aov[:-1]["sum_sq"] - (aov[:-1]["df"] * aov["mean_sq"][-1])) / (
        sum(aov["sum_sq"]) + aov["mean_sq"][-1]
    )

    cols = ["sum_sq", "df", "mean_sq", "F", "PR(>F)", "eta_sq", "omega_sq"]
    aov = aov[cols]
    return aov


def dunn_test(data, adj_method: str = "hs"):
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

    frame, groups = __array_to_df(data)

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

    result_table["Adjusted p-value"] = result_table["Adjusted p-value"].map(
        "{:.4f}".format
    )

    return result_table


def conover_iman_test(data, adj_method: str = "hs"):
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
        t = mean_rank_diff / se
        return t

    h_corr, _ = stats.kruskal(*data)
    frame, groups = __array_to_df(data)
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
    p_value_lst = list()
    for comp in comparison_lst:
        t = test_stats(frame, groups, h_corr, s2, *comp)
        p_value = 2.0 * stats.t.sf(np.abs(t), df=(total_obs - len(groups)))
        p_value_lst.append(p_value)

    hypo_test, p_adj, *_ = mt.multipletests(p_value_lst, method=adj_method)
    result_table = pd.DataFrame(
        {
            "Comparison group": comparison_lst,
            "Adjusted p-value": p_adj,
            "Reject Null?": hypo_test,
        }
    )

    result_table["Adjusted p-value"] = result_table["Adjusted p-value"].map(
        "{:.4f}".format
    )
    return result_table


if __name__ == "__main__":
    os.chdir("dna_sdr/pickles/")
    fname = "indiviual_Screen_one_phase_param.pkl"
    df = pd.read_pickle(fname)
    df.rename({"Condition": "plate_loc"}, axis=1, inplace=True)
    trig_info = pd.read_pickle("trig_info.pkl")
    trig_df = pd.merge(df, trig_info, on="plate_loc")

    # Set as future input for if specify plate/well combo
    read_lst = ["P0_A1", "P2_A5", "P2_A6"]
    trig_df.set_index("plate_loc", inplace=True)
    same_loc = trig_df.loc[read_lst]
    same_loc.reset_index(inplace=True)

    """
    Graphical repersentation of the data distribution using boxplot and catagorical scatter plot
    (swarmplot) using seaborn package.
    """
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(121)
    ax1 = sns.boxplot(
        x="plate_loc", y="rate", data=same_loc, color="#99c2a2", width=0.2
    )
    ax1 = sns.swarmplot(x="plate_loc", y="rate", data=same_loc, color="#7d0013")
    ax2 = fig.add_subplot(122)
    ax2 = sns.boxplot(
        x="plate_loc", y="plateau", data=same_loc, color="#99c2a2", width=0.2
    )
    ax2 = sns.swarmplot(x="plate_loc", y="plateau", data=same_loc, color="#7d0013")
    plt.show()

    cond = list(set([i for i in same_loc["plate_loc"]]))
    plateau_lst = list()
    rate_lst = list()
    for condition in cond:
        cond_df = same_loc[same_loc["plate_loc"] == condition]
        plateau_lst.append(cond_df["plateau"].to_numpy())
        rate_lst.append(cond_df["rate"].to_numpy())

    """
    Outlier detection using fence rule (1.5*IQR) and generalized ESD test.
    Shapiro test was also used to test the normality of the data set before and after
    the removal of outliers
    """
    for sample in range(len(plateau_lst)):
        print(cond[sample])
        sample_plateau = plateau_lst[sample]
        print("Shapiro test for the data")
        w, pvalue = stats.shapiro(sample_plateau)
        print("The W stats is {:.4f}".format(w))
        p_value_printout(pvalue)
        print("")
        fence_in, out_count = fence_rule(sample_plateau)
        print(
            "There are {} possilbe outliers based on the fence rule (1.5*IQR)".format(
                out_count
            )
        )
        if out_count != 0:
            w, pvalue = stats.shapiro(fence_in)
            print("Shapiro test after removal of potential outliers with fence rule")
            print("The W stats is {:.4f}".format(w))
            p_value_printout(pvalue)
            print("")
        print("General ESD Test for outlier")
        data_postESD = general_ESD(sample_plateau, 5, 0.05)
        print("")
        if not (np.array_equiv(data_postESD, fence_in)) or not (
            np.array_equiv(data_postESD, sample_plateau)
        ):
            w, pvalue = stats.shapiro(data_postESD)
            print("Shapiro test after removal of potential outliers with general ESD")
            print("The W stats is {:.4f}".format(w))
            p_value_printout(pvalue)
            print("")

    """
    ANOVA testing with normality plot, Shapiro test (normality test), 
    and Levene test (equal variance test) for assumpition testing.
    """
    print("")
    model = ols("rate ~ C(plate_loc)", data=same_loc).fit()
    aov = sm.stats.anova_lm(model, typ=2)
    aov_table = anova_table(aov)
    print(aov_table)

    print("")
    print("Shapiro test for normality for the ANOVA residual")
    w, pvalue = stats.shapiro(model.resid)
    print("The W stats is {:.4f}".format(w))
    p_value_printout(pvalue)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    normality_plot, stat = stats.probplot(model.resid, plot=plt, rvalue=True)
    ax.set_title("Probability plot of model's residuals", fontsize=20)
    plt.show()

    print("")

    print("Levene test for equal variance")
    levene = stats.levene(
        same_loc[same_loc["plate_loc"] == "P0_A1"]["rate"],
        same_loc[same_loc["plate_loc"] == "P2_A5"]["rate"],
        same_loc[same_loc["plate_loc"] == "P2_A6"]["rate"],
    )
    print("Test statistics is {:.4f}".format(levene[0]))
    p_value_printout(levene[-1])

    print("")
    comp = mc.MultiComparison(same_loc["rate"], same_loc["plate_loc"])
    post_hoc_res = comp.tukeyhsd()
    print(post_hoc_res.summary())
    print("")

    """
    Kruskal-Wallis H-test (non-parametric version of ANOVA) in case when the
    ANOVA assumption was not satisified. The data were also tested using 
    Kolmogorov-Smirnov test for goodness of fit to check if the samples follows
    a specific population distribution (i.e., Student's T distribution). 
    If the data follow the Student's T distribution, Conover-Iman test 
    will be used for the post-hoc analysis, else Dunn's test will be used.
    """
    print("Kruskal-Wallis H-test")
    h, p_value = stats.kruskal(*plateau_lst)
    print("the H-statistics is {:.3f}".format(h))
    p_value_printout(p_value)
