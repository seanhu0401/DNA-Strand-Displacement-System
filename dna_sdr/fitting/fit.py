"""_summary_

Returns
-------
_type_
    _description_
"""

import glob
import os

from lmfit.model import ModelResult

import pandas as pd
from lmfit import Model

from dna_sdr.fitting.model_function import (
    one_phase_association,
    lag_one_phase_association,
    sec_kin_fit,
)

from dna_sdr.experimental.data_processing import time_list_generation

time_lst = time_list_generation(60)


def one_phase_fit(data: pd.Series) -> ModelResult:
    if data.iloc[0] > 20:
        model = Model(one_phase_association)
        params = model.make_params(plateau=250, k=0.01, y_0=dict(value=20, min=0))
    else:
        model = Model(lag_one_phase_association)
        params = model.make_params(
            plateau=250, k=0.01, time_0=20, y_0=dict(value=20, min=0)
        )

    res = model.fit(data, params, time=time_lst)
    return res


def kinetic_fit(
    data: pd.Series, test: str, condition: str | None = None
) -> ModelResult:
    one_phase_result: ModelResult = one_phase_fit(data)
    plateau: float = one_phase_result.best_values["plateau"]

    if test == "Conc" and condition is not None:
        factor: float = int(condition) / 100
        y0: list[float] = [500 * factor, plateau, 0, 0]
    else:
        y0: list[float] = [500, plateau, 0, 0]

    model = Model(sec_kin_fit, independent_vars=["t", "y0"])
    params = model.make_params(k1={"value": 1e-05, "min": 1e-08, "max": 100.0})

    res = model.fit(data, params, t=time_lst, y0=y0)

    return res


def individual_fit(file: str, equation: str, labels: list[str]):
    df: pd.DataFrame = pd.read_pickle(file)
    name = file.split(".")[0]
    test = name.split("_")[2]

    cond = None

    if test == "Screen":
        trigs: str | list[str] = "_".join(name.split("_")[-2:])
        params_group_list = [trigs]
    elif test == "Conc":
        trigs = "_".join(name.split("_")[-3:-1])
        cond = name.split("_")[-1]
        if trigs.split("_")[0] == "Conc":
            trigs = name.split("_")[-2]
        params_group_list = [trigs, cond]
    elif test == "Ratio":
        trigs = [
            "_".join(name.split("_")[3:5]),
            "_".join(name.split("_")[5:7]),
        ]
        cond = "_".join(name.split("_")[-2:])
        params_group_list = [*trigs, cond]
    else:
        raise ValueError()

    params_lst: list[dict[str, str]] = []

    for count in range(len(df.columns)):
        trial = df.iloc[:, count] * 500
        param_list = params_group_list[:]
        if equation == "one_phase":
            result = one_phase_fit(trial)
            param_list.extend(
                [
                    result.best_values["plateau"],
                    result.best_values["k"],
                    result.best_values["y_0"],
                    result.rsquared,
                    result.params["plateau"].stderr,
                    result.params["k"].stderr,
                ]
            )
            if len(labels) != len(param_list):
                raise ValueError("the fit output and the equation does not match")
        elif equation == "kinetic":
            result = kinetic_fit(trial, test, cond)
            param_list.extend(
                [result.best_values["k1"], result.rsquared, result.params["k1"].stderr]
            )
            if len(labels) != len(param_list):
                raise ValueError("the fit output and the equation does not match")

        param_group_dict = dict(zip(labels, param_list))
        params_lst.append(param_group_dict)

    return params_lst


def parameter_determination(test: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The `parameter_determination` function takes in a test name and an optional parameter indicating
    whether to process individual files, and returns dataframes and lists containing parameter
    information and results for one-phase and second kinetic fits.

    Parameters
    ----------
    test : str
        The `test` parameter is a string that represents the type of test being performed. It is used to
    generate the appropriate column names for the dataframes and to determine the file query string.
    individual : bool, optional

    Returns
    -------
        The function `parameter_determination` returns a tuple containing two elements:

    """

    def column_name_generation() -> tuple[list[str], list[str]]:
        one_phase_list = [
            "plateau",
            "rate",
            "y_0",
            "r_sq_curve",
            "plateau_error",
            "rate_error",
        ]
        kinetic_list = [
            "k_rate",
            "r_sq_kin",
            "k_error",
        ]

        if test == "Screen":
            one_phase_list = ["Trig"] + one_phase_list
            kinetic_list = ["Trig"] + kinetic_list
        elif test == "Conc":
            one_phase_list = ["Trig", "Condition"] + one_phase_list
            kinetic_list = ["Trig", "Condition"] + kinetic_list
        elif test == "Ratio":
            one_phase_list = ["Trig1", "Trig2", "Condition"] + one_phase_list
            kinetic_list = ["Trig1", "Trig2", "Condition"] + kinetic_list

        return one_phase_list, kinetic_list

    def result_combination(
        param_lst: list[dict[str, str]],
        param_df_lst: list[pd.DataFrame],
    ):
        param_df = pd.DataFrame(param_lst)
        param_df_lst.append(param_df)
        return param_df_lst

    one_phase_df_lst = []
    sec_kinetic_df_lst = []

    one_phase_list, kinetic_list = column_name_generation()

    str_query = f"*_{test}_*.pkl"

    for f in glob.glob(str_query):
        print(f)

        one_phase_params_list = individual_fit(f, "one_phase", one_phase_list)
        sec_kinetic_params_list = individual_fit(f, "sec_kinetic", kinetic_list)

        one_phase_df_lst = result_combination(
            one_phase_params_list,
            one_phase_df_lst,
        )

        sec_kinetic_df_lst = result_combination(
            sec_kinetic_params_list,
            sec_kinetic_df_lst,
        )

    one_phase_params_df = pd.concat(one_phase_df_lst)
    sec_kinetic_params_df = pd.concat(sec_kinetic_df_lst)

    return one_phase_params_df, sec_kinetic_params_df


def storage(
    fit_results: pd.DataFrame,
    test: str,
    equation: str,
) -> None:
    """The `storage` function saves a DataFrame to pickle files based on the provided parameters.

    Parameters
    ----------
    fit_results : tuple[pd.DataFrame, list[dict]]
        The `fit_results` parameter is a tuple containing two elements:
    test : str
        The `test` parameter is a string that represents the type of test being performed.
    equation : str
        The "equation" parameter in the "storage" function is a string that represents the equation used
    for fitting the data. It is used to generate the file names for storing the parameter and result
    data.

    """
    new_df = fit_results
    param_fname = f"individual_{test}_{equation}_param.pkl"
    new_df.to_pickle(param_fname)


# TODO: redo the os path finding with os.path.join
def main(test: str):
    """The main function changes the current directory based on the individual parameter, calls the
    parameter_determination function with the test_type and individual parameters, changes the directory
    again, and then calls the storage function with the output and other parameters.

    Parameters
    ----------
    test : str
        The `test_type` parameter is a string that specifies the type of test being performed. It is used as an input
    to the `parameter_determination` function.

    """
    os.chdir("./dna_sdr/IO/Output/Individual/")

    output = parameter_determination(test)
    print(output)

    os.chdir("../../../..")

    os.chdir("./dna_sdr/pickles/")
    # storage(output[0], test, "one_phase")
    # storage(output[1], test, "second_kinetic")
    os.chdir("../..")


if __name__ == "__main__":
    # test_types = ["Screen", "Conc", "Ratio"]
    test_types = ["Screen"]
    for test_type in test_types:
        print(test_type)
        main(test_type)
