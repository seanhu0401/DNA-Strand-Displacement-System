"""_summary_
"""

import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

# from dna_sdr.stats.stats import general_esd


if __name__ == "__main__":
    VERSION = "v3"
    INFO = "./dna_sdr/pickles/trig_info.pkl"
    info_df: pd.DataFrame = pd.read_pickle(INFO)
    info_df = info_df.rename({"plate_loc": "Trig"}, axis=1)

    params_df: pd.DataFrame = pd.read_pickle(
        "./dna_sdr/pickles/individual_Screen_param.pkl"
    )
    params_df = params_df.replace("Screen_T1", "T1")

    analysis_df = info_df[
        (info_df["mismatch"] <= 1)
        & (info_df["toehold"] == 7)
        & (info_df["overhang"] == 5)
        & ~(info_df["name"].str.contains("T3"))
    ]
    analysis_df = analysis_df.fillna(0).astype(int, errors="ignore")
    analysis_df = analysis_df.loc[:, analysis_df.any()]

    analysis_df_with_dummies = pd.get_dummies(
        data=analysis_df,
        columns=["mismatch_type_1", "mismatch_location_1"],
        dtype=int,
        drop_first=True,
    )

    result_df = pd.merge(analysis_df_with_dummies, params_df, on="Trig")
    result_df = result_df[(result_df["r_sq_curve"] > 0.5) & (result_df["y_0"] < 300)]

    location_column_lst = [
        col for col in result_df.columns if "mismatch_location_" in col
    ]
    type_column_lst = [col for col in result_df.columns if "mismatch_type_" in col]
    full_lst = location_column_lst[1:] + type_column_lst
    trigs = list(set(result_df["Trig"]))

    T1_df = result_df[result_df["Trig"] == "T1"]
    avg_T1_plateau = np.mean(T1_df.plateau)
    avg_T1_rate = np.mean(T1_df.rate)

    result_df["plateau"] = result_df["plateau"] / avg_T1_plateau
    result_df["rate"] = result_df["rate"] / avg_T1_rate

    # for trigger in trigs:
    #     trig_series = result_df[result_df["Trig"] == trigger].reset_index()
    #     print(trig_series["plateau"])
    #     print(general_esd(trig_series["plateau"]))

    parameter = "plateau"

    location_model = " + ".join(location_column_lst)
    type_model = " + ".join(type_column_lst)
    half_type_model = f"{parameter} ~ {type_model}"
    half_location_model = f"{parameter} ~ {location_model}"
    full_model = f"{parameter} ~ {location_model} + {type_model}"

    olsr_model = smf.ols(formula=full_model, data=result_df)
    olsr_model_results = olsr_model.fit()
    summary1 = olsr_model_results.summary()
    cov_params = olsr_model_results.cov_params()

    original_name = list(cov_params.columns)
    new_name = []
    for name in original_name:
        new_name.append(name.replace("_1_", "_"))

    cov_params.index = new_name
    cov_params.columns = new_name

    summary = olsr_model_results.summary2()
    reg_result = pd.DataFrame(summary.tables[1])
    print(reg_result)

    # reg_result.to_pickle(f"./dna_sdr/pickles/dummy_reg_{parameter}_{VERSION}.pkl")
    # cov_params.to_pickle(f"./dna_sdr/pickles/dummy_reg_{parameter}_cov_{VERSION}.pkl")
    # regression_result.to_csv(f"./dna_sdr/dummy_reg_{parameter}.csv")
