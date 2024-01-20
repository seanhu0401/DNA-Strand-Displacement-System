"""
Module Docstring
"""

import pandas as pd


def fit_result_combine(
        dataframe_1: pd.DataFrame, dataframe_2: pd.DataFrame, test:str
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
                dataframe.drop(["Condition", "index", "Trig1", "Trig2"], axis=1, inplace=True)
            else:
                dataframe.drop("index", axis=1, inplace=True)
        elif test in ("Screen", "Conc"):
            if cond_drop:
                dataframe.drop(["Condition", "index"], axis=1, inplace=True)
            else:
                dataframe.drop("index", axis=1, inplace=True)
        else:
            raise ValueError()
        return dataframe

    dataframe_1 = dataframe_processing(dataframe_1)
    dataframe_2 = dataframe_processing(dataframe_2, True)
    final_df = pd.merge(dataframe_1, dataframe_2, left_index=True, right_index=True)
    if test in ("Screen", "Conc"):
        final_df.rename({"Condition": "plate_loc"}, axis=1, inplace=True)
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

    # trig_pickle = "./dna_sdr/pickles/trig_param.pkl"  # Fitting result
    # conc_pickle = "./dna_sdr/pickles/concentration.pkl"  # Nupack analysis
    # trig_info_pickle = "./dna_sdr/pickles/trig_info.pkl"  # Trigger design information

    # trig_cf_df = pd.read_pickle(trig_pickle)
    # trig_info_df = pd.read_pickle(trig_info_pickle)
    # conc_df = pd.read_pickle(conc_pickle)

    # # Reordering the columns for the dataframe
    # cols = trig_cf_df.columns.tolist()
    # cols.insert(0, "plate_loc")
    # cols.pop(6)
    # trig_cf_df = trig_cf_df[cols]

    # trig_info_cf_df = pd.merge(trig_info_df, trig_cf_df, on="plate_loc")

    # # Recover the trigger name from the nupack analysis
    # name_lst = list()
    # for index in conc_df.index:
    #     name = index.split("+")[0][1:]
    #     name_lst.append(name)

    # conc_df["name"] = name_lst
    # conc_df.reset_index(inplace=True)
    # merged_df = pd.merge(trig_info_cf_df, conc_df, on="name")
    # merged_df.drop("index", axis=1, inplace=True)
    # merged_df.rename(columns={"Conc (nM)": "Nupack"}, inplace=True)
    # merged_df["percent diff"] = (
    #     (merged_df["plateau"] - merged_df["Nupack"]) / merged_df["Nupack"] * 100
    # )

    # print(merged_df)
    # merged_df.to_csv("merged_df.csv")
