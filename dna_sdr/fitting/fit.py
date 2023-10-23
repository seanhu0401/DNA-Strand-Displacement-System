"""
xxx
"""

from math import e
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
        mean_df = df.filter(regex="mean")
        std_df = df.filter(regex="std")
        filtered_mean = mean_df.filter(regex=condition) * 500
        filtered_std = std_df.filter(regex=condition) * 500
        mean = filtered_mean.squeeze()
        std = filtered_std.squeeze()
        result = model.fit(mean, params, time=time_lst, weight=1 / std)  # type: ignore

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
        mean_df = df.filter(regex="mean")
        std_df = df.filter(regex="std")
        filtered_mean = mean_df.filter(regex=condition) * 500
        filtered_std = std_df.filter(regex=condition) * 500
        mean = filtered_mean.squeeze()
        std = filtered_std.squeeze()
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

    params_list = []
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
    else:
        raise ValueError()

    params_list = []
    full_result_dict = {}

    for condition in conditions:
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
        else:
            raise ValueError()

        param_group_dict = dict(zip(labels, params_group_list))
        params_list.append(param_group_dict)

        if not bool(full_result_dict):
            full_result_dict = result_dict
        else:
            full_result_dict.update(result_dict)

    return params_list, full_result_dict


def parameter_determination(test, individual=False):
    """
    xxx

    Parameters
    ----------
    file :

    equation :

    lable :

    Returns
    -------

    """
    one_phase_df_list = list()
    one_phase_result_dict = dict()

    sec_kinetic_df_list = list()
    sec_kinetic_result_dict = dict()
    if not individual:
        str_query = f"*_{test}_*_summerized.pkl"
        if test == "Ratio":
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

        elif test == "Conc" or test == "Screen":
            one_phase_list = [
                "Plate Number",
                "Trigger type",
                "plateau",
                "rate",
                "r_sq_curve",
                "plateau_error",
                "rate_error",
            ]

            kinetic_list = [
                "Plate Number",
                "Trigger type",
                "k_rate",
                "r_sq_kin",
                "k_error",
            ]

    else:
        str_query = f"*_{test}_*.pkl"
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

    for f in glob.glob(str_query):
        print(f)
        if not individual:
            one_phase_params_list, one_phase_results = overall_fit(
                f, test, "one_phase", one_phase_list
            )
            sec_kinetic_params_list, sec_kinetic_results = overall_fit(
                f, test, "sec_kinetic", kinetic_list
            )

        else:
            one_phase_params_list, one_phase_results = ind_fit(
                f, "one_phase", one_phase_list
            )
            sec_kinetic_params_list, sec_kinetic_results = ind_fit(
                f, "sec_kinetic", kinetic_list
            )

        if not bool(one_phase_result_dict):
            one_phase_result_dict = one_phase_results
        else:
            one_phase_result_dict.update(one_phase_results)
        one_phase_df = pd.DataFrame(one_phase_params_list)
        one_phase_df_list.append(one_phase_df)

        if not bool(sec_kinetic_result_dict):
            sec_kinetic_result_dict = sec_kinetic_results
        else:
            sec_kinetic_result_dict.update(sec_kinetic_results)
        sec_kinetic_df = pd.DataFrame(sec_kinetic_params_list)
        sec_kinetic_df_list.append(sec_kinetic_df)

    one_phase_params_df = pd.concat(one_phase_df_list)
    sec_kinetic_params_df = pd.concat(sec_kinetic_df_list)

    output = (
        one_phase_params_df,
        sec_kinetic_params_df,
        one_phase_result_dict,
        sec_kinetic_result_dict,
    )

    return output


if __name__ == "__main__":
    # os.chdir("./dna_sdr/pickles/")
    os.chdir("./dna_sdr/IO/Output/Pickles/")
    # df = pd.read_pickle("indiviual_Screen_one_phase_param.pkl")
    # df = pd.read_pickle("Ratio_one_phase_param.pkl")
    # df = pd.read_pickle("Screen_one_phase_param.pkl")
    # df = pd.read_pickle("4WJ_HEX_Screen_P0_A01_A04_summerized.pkl")
    # print(df)

    #######
    # individual = True
    # if individual:
    #     os.chdir("./dna_sdr/IO/Output/Pickles/Individual/")
    # else:
    #     os.chdir("./dna_sdr/IO/Output/Pickles/")

    # test_type = "Screen"
    # output = parameter_determination(test_type, individual=True)

    # if individual:
    #     os.chdir("../../../../..")
    # else:
    #     os.chdir("../../../..")

    # os.chdir("./dna_sdr/pickles/")

    # if individual:
    #     if os.path.exists("indiviual_{}_one_phase_param.pkl".format(test_type)):
    #         one_phase_params_df = pd.read_pickle(
    #             "indiviual_{}_one_phase_param.pkl".format(test_type)
    #         )
    #         if not one_phase_params_df.equals(output[0]):
    #             output[0].to_pickle(
    #                 "indiviual_{}_one_phase_param.pkl".format(test_type)
    #             )
    #             with open(
    #                 "indiviual_{}_one_phase_result.pkl".format(test_type), "wb"
    #             ) as fp:
    #                 pickle.dump(output[3], fp)
    #     else:
    #         output[0].to_pickle("indiviual_{}_one_phase_param.pkl".format(test_type))
    #         with open(
    #             "indiviual_{}_one_phase_result.pkl".format(test_type), "wb"
    #         ) as fp:
    #             pickle.dump(output[3], fp)

    #     if os.path.exists("indiviual_{}_first_kinetic_param.pkl".format(test_type)):
    #         first_kine_params_df = pd.read_pickle(
    #             "indiviual_{}_first_kinetic_param.pkl".format(test_type)
    #         )
    #         if not first_kine_params_df.equals(output[1]):
    #             output[1].to_pickle(
    #                 "indiviual_{}_first_kinetic_param.pkl".format(test_type)
    #             )
    #             with open(
    #                 "indiviual_{}_first_kinetic_result.pkl".format(test_type), "wb"
    #             ) as fp:
    #                 pickle.dump(output[4], fp)

    #     else:
    #         output[1].to_pickle(
    #             "indiviual_{}_first_kinetic_param.pkl".format(test_type)
    #         )
    #         with open(
    #             "indiviual_{}_first_kinetic_result.pkl".format(test_type), "wb"
    #         ) as fp:
    #             pickle.dump(output[4], fp)

    #     if os.path.exists("indiviual_{}_second_kinetic_param.pkl".format(test_type)):
    #         sec_kine_params_df = pd.read_pickle(
    #             "indiviual_{}_second_kinetic_param.pkl".format(test_type)
    #         )
    #         if not sec_kine_params_df.equals(output[2]):
    #             output[2].to_pickle(
    #                 "indiviual_{}_second_kinetic_param.pkl".format(test_type)
    #             )
    #             with open(
    #                 "indiviual_{}_second_kinetic_result.pkl".format(test_type), "wb"
    #             ) as fp:
    #                 pickle.dump(output[5], fp)

    #     else:
    #         output[2].to_pickle(
    #             "indiviual_{}_second_kinetic_param.pkl".format(test_type)
    #         )
    #         with open(
    #             "indiviual_{}_second_kinetic_result.pkl".format(test_type), "wb"
    #         ) as fp:
    #             pickle.dump(output[5], fp)

    # else:
    #     if os.path.exists("{}_one_phase_param.pkl".format(test_type)):
    #         one_phase_params_df = pd.read_pickle(
    #             "{}_one_phase_param.pkl".format(test_type)
    #         )
    #         if not one_phase_params_df.equals(output[0]):
    #             output[0].to_pickle("{}_one_phase_param.pkl".format(test_type))
    #             with open("{}_one_phase_result.pkl".format(test_type), "wb") as fp:
    #                 pickle.dump(output[3], fp)
    #     else:
    #         output[0].to_pickle("{}_one_phase_param.pkl".format(test_type))
    #         with open("{}_one_phase_result.pkl".format(test_type), "wb") as fp:
    #             pickle.dump(output[3], fp)

    #     if os.path.exists("{}_first_kinetic_param.pkl".format(test_type)):
    #         first_kine_params_df = pd.read_pickle(
    #             "{}_first_kinetic_param.pkl".format(test_type)
    #         )
    #         if not first_kine_params_df.equals(output[1]):
    #             output[2].to_pickle("{}_first_kinetic_param.pkl".format(test_type))
    #             with open("{}_first_kinetic_result.pkl".format(test_type), "wb") as fp:
    #                 pickle.dump(output[4], fp)

    #     else:
    #         output[2].to_pickle("{}_first_kinetic_param.pkl".format(test_type))
    #         with open("{}_first_kinetic_result.pkl".format(test_type), "wb") as fp:
    #             pickle.dump(output[4], fp)

    #     if os.path.exists("{}_second_kinetic_param.pkl".format(test_type)):
    #         sec_kine_params_df = pd.read_pickle(
    #             "{}_second_kinetic_param.pkl".format(test_type)
    #         )
    #         if not sec_kine_params_df.equals(output[1]):
    #             output[2].to_pickle("{}_second_kinetic_param.pkl".format(test_type))
    #             with open("{}_second_kinetic_result.pkl".format(test_type), "wb") as fp:
    #                 pickle.dump(output[5], fp)

    #     else:
    #         output[2].to_pickle("{}_second_kinetic_param.pkl".format(test_type))
    #         with open("{}_second_kinetic_result.pkl".format(test_type), "wb") as fp:
    #             pickle.dump(output[5], fp)
