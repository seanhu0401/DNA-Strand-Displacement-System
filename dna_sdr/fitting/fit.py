"""
xxx
"""
import glob
import os
import pickle

import numpy as np
import pandas as pd
from lmfit import Model
from lmfit.lineshapes import logistic
from scipy.integrate import solve_ivp

from dna_sdr.experimental.data_processing import time_list_generation

time_lst = time_list_generation(60)


def one_phase_association(time: float, plateau: float, k: float, y_0: float) -> float:
    """The function calculates the value of a one-phase association reaction over time.
    """
    return y_0 + (plateau - y_0) * (1 - np.exp(-k * time))


def lag_one_phase_association(
        time: float, plateau: float, k: float, time_0: float, y_0: float
) -> float:
    """The function calculates the value of a one-phase association reaction over time.

    Parameters
    ----------
    time : float
        The time parameter represents the time at which the association is being calculated.
    plateau : float
        The plateau is the maximum value that the function will approach as time goes to infinity.
    k : float
        The parameter "k" represents the rate constant in the one-phase association equation. It determines
    how quickly the association reaction occurs. A higher value of "k" indicates a faster association
    rate.
    time_0 : float
        The value where the reaction start after the initial lag.
    y_0 : float
        The initial value during the initial lag time.

    Returns
    -------
        y_0 if the time is less than time_0 else result of the equation
        `y_0 + (plateau - y_0) * (1 - np.exp(-k * time))`.

    """
    step = logistic(time, center=time_0, sigma=0.1)
    return (
            y_0 * (1 - step)
            + (y_0 + (plateau - y_0) * (1 - np.exp(-k * (time - time_0)))) * step
    )


def second_kinetic(t: list[float], y0: list[float], k: float) -> list[float]:
    """The function `second_kinetic` defines a simplified system of ordinary differential equations (ODEs)
    for a second order irreversible reaction.

    Parameters
    ----------
    t : list[float]
        The parameter `t` is a list of time values at which you want to evaluate the ODEs. It represents
    the time points at which you want to calculate the values of the variables in the system.
    y0 : list[float]
        The parameter `y0` is a list of initial values for the variables in the system. In this case, the
    variables are `I`, `QT`, `IT`, and `Q`. So `y0` should be a list of four floats representing the
    initial values of these variables.
    k : float
        The parameter `k` represents the rate constant for the second order irreversible reaction. It
    determines the rate at which the reaction occurs.

    Returns
    -------
        a list of floats, which represents the system of ordinary differential equations (ODEs) for the
    given second order irreversible reaction.

    """

    I, QT, IT, Q = y0
    k1 = k

    # the model equations
    dIdt = -k1 * I * QT
    dQTdt = -k1 * I * QT
    dITdt = k1 * I * QT
    dQdt = k1 * I * QT
    ode = [dIdt, dQTdt, dITdt, dQdt]
    return ode


def sec_kin_fit(t: list[float], y0: list[float], k1: float) -> np.ndarray:
    """The function `sec_kin_fit` solves a second-order kinetic equation using the `solve_ivp` function and
    returns the third element of the solution.

    Parameters
    ----------
    t : list[float]
        A list of time values at which the solution is evaluated.
    y0 : list[float]
        The parameter `y0` represents the initial conditions for the system of differential equations. In
    this case, it is a list of initial values for the dependent variables in the system. The length of
    `y0` should match the number of equations in the system.
    k1 : float
        The parameter `k1` represents the rate constant for the second-order kinetic reaction. It
    determines the rate at which the reaction proceeds.

    Returns
    -------
        The function `sec_kin_fit` returns an `np.ndarray` which represents the third element of the
    solution of the second order kinetic equation.

    """
    x = solve_ivp(
        second_kinetic,
        (0, max(t) + 10),
        y0,
        t_eval=t,
        args=(k1,),
        rtol=1e-9,
        method="RK45",
    )
    return x.y[2]


def individual_fit(file: str, equation: str, labels: list[str]):
    """The `ind_fit` function reads a DataFrame from a pickle file, performs either one-phase or second
    kinetic fitting on the data, and returns a list of parameter dictionaries and a dictionary of
    fitting results.

    Parameters
    ----------
    file : str
        The `file` parameter is a string that represents the file name or path of the pickle file that
    contains the data to be processed.
    equation : str
        _The `equation` parameter is a string that specifies the type of equation to use for fitting. It can
    have two possible values: "one_phase" or "sec_kinetic".
    labels : list[str]
        The `labels` parameter is a list of strings that represents the labels or names for the different
    parameter groups. These labels will be used to create a dictionary where each label is associated
    with a list of parameter values.

    Returns
    -------

    """

    def _one_phase_fitting(data: pd.Series, param_group_list: list[str | int]):
        if data.iloc[0] > 0:
            model = Model(one_phase_association)
            params = model.make_params(plateau=250, k=0.01, y_0=0)
        else:
            model = Model(lag_one_phase_association)
            params = model.make_params(plateau=250, k=0.01, time_0=100, y_0=0)

        res = model.fit(data, params, time=time_lst)

        param_group_list.extend(
            [
                res.best_values["plateau"],
                res.best_values["k"],
                res.rsquared,
                res.params["plateau"].stderr,
                res.params["k"].stderr,
            ]
        )
        print(param_group_list)

        if len(labels) != len(param_group_list):
            raise ValueError("the fit output and the equation does not match")

        return res, param_group_list

    def _kinetic_fitting(data: pd.Series, param_group_list: list[str | int], condition: str | None = None):
        if test == "Conc":
            factor = int(condition) / 100
            y0 = [500 * factor, 500, 0, 0]
        else:
            y0 = [500, 500, 0, 0]

        model = Model(sec_kin_fit, independent_vars=["t", "y0"])
        params = model.make_params(k1={"value": 1e-05, "min": 1e-08, "max": 100.0})

        res = model.fit(data, params, t=time_lst, y0=y0)

        param_group_list.extend(
            [res.best_values["k1"], res.rsquared, res.params["k1"].stderr]
        )

        if len(labels) != len(param_group_list):
            raise ValueError("the fit output and the equation does not match")

        return res, param_group_list

    df: pd.DataFrame = pd.read_pickle(file)
    name = file.split(".")[0]
    test = name.split("_")[2]

    if test == "Screen":
        trigs = "_".join(name.split("_")[-2:])
        params_group_list = [trigs]
    elif test == "Conc":
        trigs = "_".join(name.split("_")[-3:-1])
        cond = name.split("_")[-1]
        if trigs.split("_")[0] == "Conc":
            trigs = name.split("_")[-2]
        params_group_list = [trigs, cond]
    elif test == "Ratio":
        trigs = ["_".join(name.split("_")[3:5]), "_".join(name.split("_")[5:7])]
        cond = "_".join(name.split("_")[-2:])
        params_group_list = [*trigs, cond]
    else:
        raise ValueError()

    full_result_dict = {}
    params_list: list[dict[str, str]] = []

    for count in range(len(df.columns)):
        trial = df.iloc[:, count] * 500
        param_list = params_group_list[:]
        if equation == "one_phase":
            if test in {"Conc", "Ratio"}:
                result, param_list = _one_phase_fitting(
                    trial, param_list
                )
            else:
                result, param_list = _one_phase_fitting(
                    trial, param_list,
                )
        elif equation == "sec_kinetic":
            if test in "Conc":
                result, param_list = _kinetic_fitting(
                    trial, param_list, cond
                )
            else:
                result, param_list = _kinetic_fitting(
                    trial, param_list
                )
        else:
            raise ValueError()

        result_dict = dict({f"{count}": result})
        param_group_dict = dict(zip(labels, param_list))
        params_list.append(param_group_dict)

        if not bool(full_result_dict):
            full_result_dict = result_dict
        else:
            full_result_dict.update(result_dict)

    return params_list, full_result_dict


def parameter_determination(
        test: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, dict]], list[dict[str, dict]]]:
    """The `parameter_determination` function takes in a test name and an optional parameter indicating
    whether to process individual files, and returns dataframes and lists containing parameter
    information and results for one-phase and second kinetic fits.

    Parameters
    ----------
    test : str
        The `test` parameter is a string that represents the type of test being performed. It is used to
    generate the appropriate column names for the dataframes and to determine the file query string.
    individual : bool, optional
        The `individual` parameter is a boolean flag that determines whether the function should perform
    individual fitting or overall fitting. If `individual` is set to `True`, the function will perform
    individual fitting. If `individual` is set to `False`, the function will perform overall fitting.
    The default value for

    Returns
    -------
        The function `parameter_determination` returns a tuple containing four elements:

    """

    def column_name_generation() -> tuple[list[str], list[str]]:
        one_phase_list = [
            "plateau",
            "rate",
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
            file: str,
            param_lst: list[dict[str, str]],
            result: dict,
            result_lst: list[dict],
            param_df_lst: list[pd.DataFrame],
    ):
        test = file.split(".")[0].split("Conc_")[-1]
        result_lst.append({f"{test}": result})
        param_df = pd.DataFrame(param_lst)
        param_df_lst.append(param_df)
        return result_lst, param_df_lst

    one_phase_df_lst = []
    one_phase_result_lst = []

    sec_kinetic_df_lst = []
    sec_kinetic_result_lst = []
    one_phase_list, kinetic_list = column_name_generation()

    str_query = f"*_{test}_*.pkl"

    for f in glob.glob(str_query):
        print(f)

        one_phase_params_list, one_phase_results = individual_fit(
            f, "one_phase", one_phase_list
        )
        sec_kinetic_params_list, sec_kinetic_results = individual_fit(
            f, "sec_kinetic", kinetic_list
        )

        one_phase_result_lst, one_phase_df_lst = result_combination(
            f,
            one_phase_params_list,
            one_phase_results,
            one_phase_result_lst,
            one_phase_df_lst,
        )

        sec_kinetic_result_lst, sec_kinetic_df_lst = result_combination(
            f,
            sec_kinetic_params_list,
            sec_kinetic_results,
            sec_kinetic_result_lst,
            sec_kinetic_df_lst,
        )

    one_phase_params_df = pd.concat(one_phase_df_lst)
    sec_kinetic_params_df = pd.concat(sec_kinetic_df_lst)

    return (
        one_phase_params_df,
        sec_kinetic_params_df,
        one_phase_result_lst,
        sec_kinetic_result_lst,
    )


def storage(
        fit_results: tuple[pd.DataFrame, list[dict[str, dict]]],
        test: str,
        equation: str,
) -> None:
    """The `storage` function saves a DataFrame and a dictionary to pickle files based on the provided
    parameters.

    Parameters
    ----------
    fit_results : tuple[pd.DataFrame, list[dict]]
        The `fit_results` parameter is a tuple containing two elements:
    individual : bool
        The `individual` parameter is a boolean flag that indicates whether the storage should be done for
    individual results or not. If `individual` is `True`, the filenames for the stored data will include
    the prefix "individual_", otherwise they will not.
    test_type : str
        The `test_type` parameter is a string that represents the type of test being performed. It could be
    something like "regression", "classification", or any other type of analysis.
    equation : str
        The "equation" parameter in the "storage" function is a string that represents the equation used
    for fitting the data. It is used to generate the file names for storing the parameter and result
    data.

    """
    new_df, result_dict = fit_results
    param_fname = f"individual_{test}_{equation}_param.pkl"
    result_fname = f"individual_{test}_{equation}_result.pkl"
    new_df.to_pickle(param_fname)
    with open(result_fname, "wb") as fp:
        pickle.dump(result_dict, fp)


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

    output = parameter_determination(test_type)

    os.chdir("../../../..")

    os.chdir("./dna_sdr/pickles/")
    one_phase_results = (output[0], output[2])
    storage(one_phase_results, test, "one_phase")

    kinetic_results = (output[1], output[3])
    storage(kinetic_results, test, "second_kinetic")
    os.chdir("../..")


if __name__ == "__main__":
    # test_types = ["Screen", "Conc", "Ratio"]
    test_types = ["Screen"]
    for test_type in test_types:
        print(test_type)
        main(test_type)
