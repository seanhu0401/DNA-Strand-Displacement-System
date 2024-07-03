"""_summary_
"""

import numpy as np
import pandas as pd
import scipy.stats


if __name__ == "__main__":
    trig_info_df = pd.read_pickle("./dna_sdr/pickles/trig_info.pkl")
    trig_info_df = trig_info_df.rename({"plate_loc": "Trig"}, axis=1)
    # target_trig_df = (
    #     trig_info_df[
    #         (trig_info_df["name"].str.contains("T3")) & (trig_info_df["mismatch"] > 0)
    #     ]
    #     .dropna(axis=1, how="all")
    #     .fillna(0)
    #     .astype(int, errors="ignore")
    # )
    # target_trig_df = target_trig_df.replace("P3_A00", "T3")

    target_trig_df = (
        trig_info_df[(trig_info_df["toehold"] == 7) & (trig_info_df["mismatch"] == 2)]
        .dropna(axis=1, how="all")
        .fillna(0)
        .astype(int, errors="ignore")
    )

    EXPERIMENTAL = "./dna_sdr/pickles/individual_Screen_param.pkl"
    exp_df: pd.DataFrame = pd.read_pickle(EXPERIMENTAL)
    exp_df = exp_df.replace({"Screen_T1": "T1", "Screen_T3": "T3"})
    target_df = pd.merge(target_trig_df, exp_df, on="Trig")
    curve_df = target_df[target_df["r_sq_curve"] > 0.5]

    T3_df = exp_df[exp_df["Trig"] == "T3"]
    avg_T3_plateau = np.mean(T3_df["plateau"])
    avg_T3_rate = np.mean(T3_df["rate"])

    curve_df["plateau"] = curve_df["plateau"] / avg_T3_plateau
    curve_df["rate"] = curve_df["rate"] / avg_T3_rate

    curve_group_df = pd.DataFrame(
        [
            curve_df.groupby(["seq"])["plateau"].mean(),
            curve_df.groupby(["seq"])["plateau"].median(),
            curve_df.groupby(["seq"])["plateau"].std(),
            curve_df.groupby(["seq"])["rate"].mean(),
            curve_df.groupby(["seq"])["rate"].median(),
            curve_df.groupby(["seq"])["rate"].std(),
        ],
        index=[
            "plateau_mean",
            "platuea_median",
            "plateau_std",
            "rate_mean",
            "rate_median",
            "rate_std",
        ],
    ).T.reset_index()

    mismatch_loc_lst = []
    mismatch_type_lst = []

    for i in range(len(target_trig_df)):
        mismatch_loc_lst.append(
            [
                "mismatch_location_" + str(loc)
                for loc in target_trig_df["mismatch_loc"].iloc[i]
            ]
        )

        mismatch_type_lst.append(
            [
                "mismatch_type_" + str(types)
                for types in target_trig_df["mismatch_type"].iloc[i]
            ]
        )

    seq_lst = list(target_trig_df["seq"])

    PARMAS = "plateau"
    VERSION = "v3"
    parameter_df = pd.read_pickle(
        f"./dna_sdr/pickles/dummy_reg_{PARMAS}_{VERSION}.pkl"
    ).reset_index()
    parameter_df["index"] = parameter_df["index"].str.replace("_1_", "_")
    parameter_df = parameter_df.set_index("index")

    p_values_array = parameter_df.iloc[:, 3]

    cov_df = pd.read_pickle(f"./dna_sdr/pickles/dummy_reg_{PARMAS}_cov_{VERSION}.pkl")

    intercept = parameter_df.loc["Intercept"].iloc[0]

    mismatch_df = pd.DataFrame(
        [seq_lst, mismatch_loc_lst, mismatch_type_lst],
        index=["seq", "mismatch_location", "mismatch_type"],
    )

    t_crit = scipy.stats.t.ppf(q=1 - 0.05 / 2, df=595)
    pred_lst = []
    se_lst = []

    for _, row in mismatch_df.items():
        mismatch = []
        mismatch.extend(row.iloc[1])
        mismatch.extend(row.iloc[2])
        se = np.sqrt(cov_df[mismatch].loc[mismatch].sum().sum())

        # se = parameter_df.loc["Intercept"].iloc[1]

        # if (p_values_array[row.iloc[1]].iloc[0] < 0.05) & (
        #     p_values_array[row.iloc[2]].iloc[0] < 0.05
        # ):
        #     se = np.sqrt(cov_df[mismatch].loc[mismatch].sum().sum())
        # elif p_values_array[row.iloc[1]].iloc[0] < 0.05:
        #     se = np.sqrt(cov_df[row.iloc[1]].loc[row.iloc[1]].sum().sum())
        # elif p_values_array[row.iloc[2]].iloc[0] < 0.05:
        #     se = np.sqrt(cov_df[row.iloc[2]].loc[row.iloc[2]].sum().sum())

        se_lst.append(se)

        mismatch_loc = parameter_df.loc[row.iloc[1]]["Coef."].sum()
        mismatch_type = parameter_df.loc[row.iloc[2]]["Coef."].sum()

        # mismatch_loc = 0
        # mismatch_type = 0

        # if p_values_array[row.iloc[1]].iloc[0] < 0.05:
        #     mismatch_loc = parameter_df.loc[row.iloc[1]]["Coef."].sum()
        # if p_values_array[row.iloc[2]].iloc[0] < 0.05:
        #     mismatch_type = parameter_df.loc[row.iloc[2]]["Coef."].sum()

        pred_value = mismatch_loc + mismatch_type + intercept
        pred_lst.append(pred_value)

    pred_df = pd.concat(
        [
            mismatch_df.T,
            pd.DataFrame({f"Prediction_{PARMAS}": pred_lst, f"SE_{PARMAS}": se_lst}),
        ],
        axis=1,
    )
    pred_df[f"0.025_{PARMAS}"] = (
        pred_df[f"Prediction_{PARMAS}"] - pred_df[f"SE_{PARMAS}"] * t_crit
    )
    pred_df[f"0.975_{PARMAS}"] = (
        pred_df[f"Prediction_{PARMAS}"] + pred_df[f"SE_{PARMAS}"] * t_crit
    )
    pred_df[f"range_{PARMAS}"] = pred_df[f"SE_{PARMAS}"] * t_crit

    res_df = target_trig_df.merge(
        pred_df,
        how="inner",
        left_on=["seq"],
        right_on=["seq"],
    )

    res_df = res_df.drop(
        labels=[
            "mismatch_location",
            "mismatch_type_x",
            "mismatch_type_y",
        ],
        axis=1,
    )

    final_df = res_df.merge(curve_group_df, on="seq")
    print(final_df[f"Prediction_{PARMAS}"])
    # final_df.to_csv(f"{PARMAS}_dummy_regression_{VERSION}.csv")
    # final_df.to_csv(f"{PARMAS}_dummy_regression_v5.csv")
    # final_df.to_pickle(f"./dna_sdr/pickles/{PARMAS}_dummy_regression_{VERSION}.pkl")
