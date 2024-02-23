"""_summary_
"""

import pandas as pd
import statsmodels.formula.api as smf


if __name__ == "__main__":
    INFO = "./dna_sdr/pickles/trig_info.pkl"
    info_df: pd.DataFrame = pd.read_pickle(INFO)
    info_df = info_df.rename({"plate_loc": "Trig"}, axis=1)
    print(info_df["Trig"])

    analysis_df = info_df[
        (info_df["mismatch"] <= 1)
        & (info_df["toehold"] == 7)
        & (info_df["overhang"] == 5)
    ]
    analysis_df = analysis_df.fillna(0).astype(int, errors="ignore")
    analysis_df = analysis_df.loc[:, analysis_df.any()]

    analysis_df_with_dummies = pd.get_dummies(
        data=analysis_df,
        columns=["mismatch type 1", "mismatch location 1"],
        dtype=int,
        drop_first=True,
    )

    location_column_lst = [
        col for col in analysis_df_with_dummies.columns if "mismatch location" in col
    ]
    type_column_lst = [
        col for col in analysis_df_with_dummies.columns if "mismatch type" in col
    ]
