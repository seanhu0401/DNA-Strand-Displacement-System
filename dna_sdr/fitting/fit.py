import pandas as pd
import numpy as np
from lmfit import Model
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os, glob, pickle


def one_phase_association(time, plateau, k):
    # pseudo-first order association kinetics
    return plateau * (1 - np.exp(-k * time))


def kinetic(t, y0, k1):
    # simplified system of ODE for SDR system
    I, QT, IT, Q = y0
    k1 = k1

    # the model equations
    dIdt = -k1 * I * QT
    dQTdt = -k1 * I * QT
    dITdt = k1 * I * QT
    dQdt = k1 * I * QT
    ode = [dIdt, dQTdt, dITdt, dQdt]
    return ode


def kin_fit(t, y0, k1):
    x = solve_ivp(
        kinetic,
        (1, max(t) + 10),
        y0,
        t_eval=t,
        args=(k1,),
        rtol=1e-9,
        method="LSODA",
    )
    return x.y[2]


def fit(file, test, equation: str, labels: list):
    df = pd.read_pickle(file)
    plate_num = file.split("_")[3]
    time_df = df["time (min)"]
    mean_df = df.filter(regex="mean")
    std_df = df.filter(regex="std")
    y0 = [500, 500, 0, 0]
    if test == "screen":
        conditions = [i.split("_")[0] for i in mean_df.columns]
    elif test == "Conc":
        if file.split("_")[4] == "A00":
            conditions = ["_".join(i.split("_")[0:-1]) for i in mean_df.columns]
        else:
            conditions = ["_".join(i.split("_")[1:-1]) for i in mean_df.columns]
            conditions[0] = "T1"

    params_list = list()
    full_result_dict = dict()

    for condition in conditions:
        params_group_list = [plate_num, condition]
        mean = mean_df.filter(regex=condition) * 500
        std = std_df.filter(regex=condition) * 500
        mean = mean.squeeze()
        std = std.squeeze()

        if equation == "one_phase":
            model = Model(one_phase_association)
            params = model.make_params(plateau=250, k=0.01)
            result = model.fit(mean, params, time=time_df, weights=1 / std)
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

        elif equation == "kinetic":
            model = Model(kin_fit, independent_vars=["t", "y0"])
            params = model.make_params(k1=dict(value=1e-5, min=10e-9, max=10e1))
            result = model.fit(mean, params, t=time_df, y0=y0, weights=1 / std)
            params_group_list.extend(
                [result.params["k1"].value, result.rsquared, result.params["k1"].stderr]
            )
            if len(labels) != len(params_group_list):
                raise ValueError("the fit output and the equation does not match")

        result_dict = dict({"{}_{}".format(plate_num, condition): result})

        param_group_dict = dict(zip(labels, params_group_list))
        params_list.append(param_group_dict)
        if not bool(full_result_dict):
            full_result_dict = result_dict
        else:
            full_result_dict.update(result_dict)

    return params_list, full_result_dict


def parameter_determination(test):
    one_phase_df_list = list()
    one_phase_result_dict = dict()
    one_phase_list = [
        "Plate Number",
        "Trigger type",
        "plateau",
        "rate",
        "r_sq_curve",
        "plateau_error",
        "rate_error",
    ]

    kinetic_df_list = list()
    kinetic_result_dict = dict()
    kinetic_list = [
        "Plate Number",
        "Trigger type",
        "k_rate",
        "r_sq_kin",
        "k_error",
    ]

    for f in glob.glob("*" + test + "*" + "summerized.pkl"):
        print(f)
        one_phase_params_list, one_phase_results = fit(
            f, test, "one_phase", one_phase_list
        )
        if not bool(one_phase_result_dict):
            one_phase_result_dict = one_phase_results
        else:
            one_phase_result_dict.update(one_phase_results)
        one_phase_df = pd.DataFrame(one_phase_params_list)
        one_phase_df_list.append(one_phase_df)

        kinetic_params_list, kinetic_results = fit(f, test, "kinetic", kinetic_list)
        if not bool(kinetic_result_dict):
            kinetic_result_dict = kinetic_results
        else:
            kinetic_result_dict.update(kinetic_results)
        kinetic_df = pd.DataFrame(kinetic_params_list)
        kinetic_df_list.append(kinetic_df)

    one_phase_params_df = pd.concat(one_phase_df_list)
    kinetic_params_df = pd.concat(kinetic_df_list)

    output = (
        one_phase_params_df,
        kinetic_params_df,
        one_phase_result_dict,
        kinetic_result_dict,
    )

    return output


if __name__ == "__main__":
    os.chdir("./dna_sdr/IO/Output/Pickles")
    test_type = "Conc"

    output = parameter_determination(test_type)

    os.chdir("../../../..")
    os.chdir("./dna_sdr/pickles/")

    if os.path.exists("{}_one_phase_param.pkl".format(test_type)):
        one_phase_params_df = pd.read_pickle("{}_one_phase_param.pkl".format(test_type))
        if not one_phase_params_df.equals(output[0]):
            output[0].to_pickle("{}_one_phase_param.pkl".format(test_type))
            with open("{}_one_phase_result.pkl".format(test_type), "wb") as fp:
                pickle.dump(output[2], fp)
    else:
        output[0].to_pickle("{}_one_phase_param.pkl".format(test_type))
        with open("{}_one_phase_result.pkl".format(test_type), "wb") as fp:
            pickle.dump(output[2], fp)

    if os.path.exists("{}_kinetic_param.pkl".format(test_type)):
        kine_params_df = pd.read_pickle("{}_kinetic_param.pkl".format(test_type))
        if not kine_params_df.equals(output[1]):
            output[1].to_pickle("{}_kinetic_param.pkl".format(test_type))
            with open("{}_kinetic_result.pkl".format(test_type), "wb") as fp:
                pickle.dump(output[3], fp)

    else:
        output[1].to_pickle("{}_kinetic_param.pkl".format(test_type))
        with open("{}_kinetic_result.pkl".format(test_type), "wb") as fp:
            pickle.dump(output[3], fp)
