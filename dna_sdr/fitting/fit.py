"""
xxx
"""

import os
import glob
import pickle
import pandas as pd
import numpy as np
from lmfit import Model
from scipy.integrate import solve_ivp
from dna_sdr.experimental.data_processing import time_list_generation


def one_phase_association(time: float, plateau: float, k: float) -> float:
    """
    Define pseudo-first order association kinetics

    Parameters
    ----------
    time : float

    plateau : float

    k : float

    Returns
    -------

    """
    return plateau * (1 - np.exp(-k * time))


def second_kinetic(t: list[float], y0: list[float], k: float) -> list[float]:
    """
    Define simplified system of ODE for SDR system

    Second order irreversible reaction:
    I + QT -> IT + Q

    Parameters
    ----------
    t : list[float]

    y0 : list[float]

    k : float

    Returns
    -------
    ode : list[float]


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
    """
    Fitting second order kinetic using IVP solver

    Parameters
    ----------
    t : list[float]

    y0 : list[float]

    k1 : float

    Returns
    -------


    """
    x = solve_ivp(
        second_kinetic,
        (1, max(t) + 10),
        y0,
        t_eval=t,
        args=(k1,),
        rtol=1e-9,
        method="LSODA",
    )
    return x.y[2]


def _one_phase_fitting(
    df: pd.DataFrame | pd.Series,
    condition: str,
    labels: list[str],
    individual_fit: bool = False,
):
    time_lst = time_list_generation(60)
    trial = df
    params_group_list = [condition]
    model = Model(one_phase_association)
    params = model.make_params(plateau=250, k=0.01)

    if not individual_fit:
        mean = (df.filter(regex="mean").filter(regex=condition) * 500).squeeze()
        std = (df.filter(regex="std").filter(regex=condition) * 500).squeeze()
        result = model.fit(mean, params, t=time_lst, weight=1 / std)  # type: ignore

    else:
        result = model.fit(trial, params, time=time_lst)

    params_group_list.extend(
        [
            result.values["plateau"],
            result.values["k"],
            result.rsquared,
            result.params["plateau"].stderr,
            result.params["k"].stderr,
        ]
    )

    if len(labels) != len(params_group_list):
        raise ValueError("the fit output and the equation does not match")

    return result, params_group_list


def _kinetic_fitting(
    df: pd.DataFrame | pd.Series,
    condition: str,
    labels: list[str],
    individual_fit: bool = False,
):
    y0 = [500, 500, 0, 0]
    time_lst = time_list_generation(60)
    trial = df
    params_group_list = [condition]
    model = Model(sec_kin_fit, independent_vars=["t", "y0"])
    params = model.make_params(k1={"value": 1e-05, "min": 1e-08, "max": 100.0})

    if not individual_fit:
        mean = (df.filter(regex="mean").filter(regex=condition) * 500).squeeze()
        std = (df.filter(regex="std").filter(regex=condition) * 500).squeeze()
        result = model.fit(mean, params, t=time_lst, y0=y0, weight=1 / std)  # type: ignore
    else:
        result = model.fit(trial, params, t=time_lst, y0=y0)

    params_group_list.extend(
        [result.params["k1"].value, result.rsquared, result.params["k1"].stderr]
    )
    if len(labels) != len(params_group_list):
        raise ValueError("the fit output and the equation does not match")

    return result, params_group_list


def ind_fit(file: str, equation: str, labels: list[str]):
    """
    xxx

    Parameters
    ----------
    file : str

    equation : str

    lable : list[str]

    Returns
    -------
    params_list :

    full_result_dict :

    """
    df: pd.DataFrame = pd.read_pickle(
        file
    )  # read the stored pickle file into dataframe
    name = file.split(".")[0]
    cond = "_".join(name.split("_")[-2:])

    params_list: list[dict[str, str]] = []
    full_result_dict = {}

    for count in range(len(df.columns)):
        trial = df.iloc[:, count] * 500
        if equation == "one_phase":
            result, params_group_list = _one_phase_fitting(trial, cond, labels, True)
        elif equation == "sec_kinetic":
            result, params_group_list = _kinetic_fitting(trial, cond, labels, True)
        else:
            raise ValueError()

        result_dict = dict({f"{count}": result})
        param_group_dict = dict(zip(labels, params_group_list))

        params_list.append(param_group_dict)
        if not bool(full_result_dict):
            full_result_dict = result_dict
        else:
            full_result_dict.update(result_dict)

    return params_list, full_result_dict


def overall_fit(file: str, test: str, equation: str, labels: list[str]):
    """
    xxx

    Parameters
    ----------
    file :

    test :

    equation :

    lable :

    Returns
    -------
    params_list :

    full_result_dict :

    """
    df: pd.DataFrame = pd.read_pickle(file)
    plate_num = file.split("_")[3]
    mean_df = df.filter(regex="mean")

    if test not in ("Ratio", "Screen", "Conc"):
        raise ValueError("Only valid test methods are 'Conc', 'Screen', and 'Ratio'")

    if test == "Screen":
        conditions = [i.split("_")[0] for i in mean_df.columns]
    elif test == "Conc":
        if file.split("_")[4] == "A00":
            conditions = ["_".join(i.split("_")[0:-1]) for i in mean_df.columns]
        else:
            conditions = ["_".join(i.split("_")[1:-1]) for i in mean_df.columns]
            conditions[0] = "T1"
    elif test == "Ratio":
        conditions = ["_".join(i.split("_")[0:-1]) for i in mean_df.columns]

    params_list: list[dict[str, str]] = []
    full_result_dict = {}

    for condition in conditions:  # type: ignore - First if statement should ensure it's not unbound
        if equation == "one_phase":
            result, params_group_list = _one_phase_fitting(df, condition, labels)
        elif equation == "sec_kinetic":
            result, params_group_list = _kinetic_fitting(df, condition, labels)
        else:
            raise ValueError()

        if test != "Ratio":
            params_group_list.insert(0, plate_num)
            result_dict = dict({f"{plate_num}_{condition}": result})
        elif test == "Ratio":
            result_dict = dict({f"{condition}": result})

        param_group_dict = dict(zip(labels, params_group_list))
        params_list.append(param_group_dict)

        if not bool(full_result_dict):
            full_result_dict = result_dict
        else:
            full_result_dict.update(result_dict)

    return params_list, full_result_dict


def parameter_determination(
    test: str, individual: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    """
    xxx

    Parameters
    ----------
    test : str

    individual : bool

    Returns
    -------

    """

    def column_name_generation(individual: bool) -> tuple[str, list[str], list[str]]:
        one_phase_list = [
            "Condition",
            "plateau",
            "rate",
            "r_sq_curve",
            "plateau_error",
            "rate_error",
        ]
        kinetic_list = [
            "Condition",
            "k_rate",
            "r_sq_kin",
            "k_error",
        ]
        if individual:
            str_query = f"*_{test}_*.pkl"
        else:
            str_query = f"*_{test}_*_summerized.pkl"
            if test in ("Conc", "Screen"):
                one_phase_list.insert(0, "Plate Number")
                kinetic_list.insert(0, "Plate Number")
            else:
                raise ValueError()

        return str_query, one_phase_list, kinetic_list

    def result_generation(file: str, eq: str, col_list: list[str], individual: bool):
        if individual:
            param_lst, result = ind_fit(file, eq, col_list)
        else:
            param_lst, result = overall_fit(file, test, eq, col_list)
        return param_lst, result

    def result_combination(
        param_lst: list[dict[str, str]],
        result: dict,
        result_dict: dict,
        param_df_lst: list[pd.DataFrame],
    ):
        if not bool(result_dict):
            result_dict = result
        else:
            result_dict.update(result)
        param_df = pd.DataFrame(param_lst)
        param_df_lst.append(param_df)

        return result_dict, param_df_lst

    one_phase_df_lst = []
    one_phase_result_dict = {}

    sec_kinetic_df_lst = []
    sec_kinetic_result_dict = {}
    str_query, one_phase_list, kinetic_list = column_name_generation(individual)

    for f in glob.glob(str_query):
        print(f)
        one_phase_params_list, one_phase_results = result_generation(
            f, "one_phase", one_phase_list, individual
        )
        sec_kinetic_params_list, sec_kinetic_results = result_generation(
            f, "sec_kinetic", kinetic_list, individual
        )

        one_phase_result_dict, one_phase_df_lst = result_combination(
            one_phase_params_list,
            one_phase_results,
            one_phase_result_dict,
            one_phase_df_lst,
        )

        sec_kinetic_result_dict, sec_kinetic_df_lst = result_combination(
            sec_kinetic_params_list,
            sec_kinetic_results,
            sec_kinetic_result_dict,
            sec_kinetic_df_lst,
        )

    one_phase_params_df = pd.concat(one_phase_df_lst)
    sec_kinetic_params_df = pd.concat(sec_kinetic_df_lst)

    return (
        one_phase_params_df,
        sec_kinetic_params_df,
        one_phase_result_dict,
        sec_kinetic_result_dict,
    )


def storage(
    fit_results: tuple[pd.DataFrame, dict],
    individual: bool,
    test_type: str,
    equation: str,
) -> None:
    """
    xxx

    Parameters
    ----------
    fit_results : tuple[pd.DataFrame, dict[str, ModelResult]]

    individual : bool

    test_type : str

    equation : str

    """
    new_df, result_dict = fit_results
    if individual:
        param_fname = f"individual_{test_type}_{equation}_param.pkl"
        result_fname = f"individual_{test_type}_{equation}_result.pkl"
    else:
        param_fname = f"{test_type}_{equation}_param.pkl"
        result_fname = f"{test_type}_{equation}_result.pkl"

    if os.path.exists(param_fname):
        current_df: pd.DataFrame = pd.read_pickle(param_fname)
        if not current_df.equals(new_df):
            new_df.to_pickle(param_fname)
            with open(result_fname, "wb") as fp:
                pickle.dump(result_dict, fp)


def main(individual: bool, test_type: str):
    """
    xxx

    Parameters
    ----------
    test : str

    individual : bool

    Returns
    -------

    """
    if individual:
        os.chdir("./dna_sdr/IO/Output/Pickles/Individual/")
    else:
        os.chdir("./dna_sdr/IO/Output/Pickles/")

    output = parameter_determination(test_type, individual)

    if individual:
        os.chdir("../../../../..")
    else:
        os.chdir("../../../..")

    os.chdir("./dna_sdr/pickles/")
    one_phase_results = (output[0], output[2])
    storage(one_phase_results, individual, test_type, "one_phase")

    kinetic_results = (output[1], output[3])
    storage(kinetic_results, individual, test_type, "second_kinetic")


if __name__ == "__main__":
    INDIVIDUAL = True
    TEST_TYPE = "Screen"
    main(INDIVIDUAL, TEST_TYPE)
