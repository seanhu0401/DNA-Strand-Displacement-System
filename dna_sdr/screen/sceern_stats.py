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


if __name__ == "__main__":

    toehold = "toehold"
    toehold_label = "Toehold (nt)"

    mis_type = "mismatch_type_1"
    mis_type_label = "Mismatch Type"

    mis_loc = "mismatch_location_1"
    mis_loc_label = "Mismatch Location"

    VERSION = "V1"
    INFO = "./dna_sdr/pickles/trig_info.pkl"
    info_df: pd.DataFrame = pd.read_pickle(INFO)
    info_df = info_df.rename({"plate_loc": "Trig"}, axis=1)
    PARMAS = "./dna_sdr/pickles/individual_Screen_param.pkl"
    params_df: pd.DataFrame = pd.read_pickle(PARMAS)
    result_df = pd.merge(info_df, params_df, on="Trig")
    result_df["mismatch_location_1"] = pd.to_numeric(
        result_df.mismatch_location_1, downcast="integer"
    )
    result_df["mismatch_type_1"] = pd.to_numeric(
        result_df.mismatch_type_1, downcast="integer"
    )

    READ = ["P0_A5", "P2_A1", "P2_A5", "P2_A10", "P2_D9", "P3_A1", "P3_A5"]
    read_df = result_df[result_df["Trig"].isin(READ)]
    read_df["rate"] = np.log(read_df["rate"])
    read_df["k_rate"] = np.log(read_df["k_rate"])
    curve_df = read_df[read_df["r_sq_curve"] >= 0.5]
    kinetic_df = read_df[read_df["r_sq_kin"] >= 0.5]
    type_2 = curve_df[curve_df[mis_type] == 2]
    type_8 = curve_df[curve_df[mis_type] == 8]

    test_param = ["k_rate", "plateau", "rate"]
    result_lst = []
    for params in test_param:
        ttest = [params]
        result = list(ttest_ind(type_2[params], type_8[params], usevar="unequal"))
        ttest.extend(result)  # type: ignore
        # ttest.append(cohens_d(type_2[params], type_8[params]))  # type: ignore
        result_lst.append(ttest)

    ttest_df = pd.DataFrame(
        result_lst,
        columns=["params", "tstats", "pvalue", "df"],
    )
    print(ttest_df)
