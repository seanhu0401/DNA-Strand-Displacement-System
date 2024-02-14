"""
Module Docstring
"""

import pandas as pd


def fit_result_combine(
    dataframe_1: pd.DataFrame, dataframe_2: pd.DataFrame, test: str
) -> pd.DataFrame:
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """

    def dataframe_processing(
        dataframe: pd.DataFrame, cond_drop: bool = False
    ) -> pd.DataFrame:
        dataframe.reset_index(inplace=True)
        if test == "Ratio":
            if cond_drop:
                dataframe.drop(
                    ["Condition", "index", "Trig1", "Trig2"], axis=1, inplace=True
                )
            else:
                dataframe.drop("index", axis=1, inplace=True)
        elif test == "Conc":
            if cond_drop:
                dataframe.drop(["Trig", "Condition", "index"], axis=1, inplace=True)
            else:
                dataframe.drop("index", axis=1, inplace=True)
        elif test == "Screen":
            if cond_drop:
                dataframe.drop(["Trig", "index"], axis=1, inplace=True)
            else:
                dataframe.drop("index", axis=1, inplace=True)
        else:
            raise ValueError()
        return dataframe

    dataframe_1 = dataframe_processing(dataframe_1)
    dataframe_2 = dataframe_processing(dataframe_2, True)
    final_df = pd.merge(dataframe_1, dataframe_2, left_index=True, right_index=True)
    if test in ("Conc", "Ratio"):
        final_df.rename({"Condition": "Concentration"}, axis=1, inplace=True)
    return final_df


if __name__ == "__main__":
    TEST = "Ratio"
    FIT1 = "one_phase"
    FIT2 = "second_kinetic"
    PARAMS1 = f"./dna_sdr/pickles/individual_{TEST}_{FIT1}_param.pkl"
    PARAMS2 = f"./dna_sdr/pickles/individual_{TEST}_{FIT2}_param.pkl"
    params_df_1: pd.DataFrame = pd.read_pickle(PARAMS1)
    params_df_2: pd.DataFrame = pd.read_pickle(PARAMS2)
    params_df = fit_result_combine(params_df_1, params_df_2, TEST)
    params_df.to_pickle(f"./dna_sdr/pickles/individual_{TEST}_param.pkl")
    # params_df.to_csv(f"./dna_sdr/pickles/individual_{TEST}_param.csv")
