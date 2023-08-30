import pandas as pd
import numpy as np
from lmfit import Model
from scipy.integrate import solve_ivp
from dna_sdr.data_process.list_generation import time_list_generation
import os, glob, pickle


def one_phase_association(time, plateau, k):
    # pseudo-first order association kinetics
    return plateau * (1 - np.exp(-k * time))


def first_kinetic(t, k1):
    return np.exp(-k1 * t)


def second_kinetic(t, y0, k1):
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
        second_kinetic,
        (1, max(t) + 10),
        y0,
        t_eval=t,
        args=(k1,),
        rtol=1e-9,
        method="LSODA",
    )
    return x.y[2]


def ind_fit(file, test, equation: str, labels: list):
    df = pd.read_pickle(file)
    name = file.split(".")[0]
    cond = "_".join(name.split("_")[-2:])
    # plate_num = name.split("_")[-2]
    # cond = name.split("_")[-1]

    time_lst = time_list_generation(60)
    y0 = [500, 500, 0, 0]

    params_list = list()
    full_result_dict = dict()

    for count in range(len(df.columns)):
        params_group_list = [cond]
        trial = df.iloc[:, count]
        if equation == "one_phase":
            model = Model(one_phase_association)
            params = model.make_params(plateau=250, k=0.01)
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

        elif equation == "first_kinetic":
            model = Model(first_kinetic)
            params = model.make_params(k1=dict(value=1e-2, min=10e-9, max=10e3))
            result = model.fit(trial, params, t=time_lst)
            params_group_list.extend(
                [result.params["k1"].value, result.rsquared, result.params["k1"].stderr]
            )
            # print(result.rsquared)
            if len(labels) != len(params_group_list):
                raise ValueError("the fit output and the equation does not match")

        elif equation == "sec_kinetic":
            model = Model(kin_fit, independent_vars=["t", "y0"])
            params = model.make_params(k1=dict(value=1e-5, min=10e-9, max=10e1))
            result = model.fit(trial, params, t=time_lst, y0=y0)
            params_group_list.extend(
                [result.params["k1"].value, result.rsquared, result.params["k1"].stderr]
            )
            print(result.rsquared)
            if len(labels) != len(params_group_list):
                raise ValueError("the fit output and the equation does not match")

        result_dict = dict({"{}".format(count): result})
        param_group_dict = dict(zip(labels, params_group_list))
        params_list.append(param_group_dict)
        if not bool(full_result_dict):
            full_result_dict = result_dict
        else:
            full_result_dict.update(result_dict)

    return params_list, full_result_dict


def overall_fit(file, test, equation: str, labels: list):
    df = pd.read_pickle(file)
    plate_num = file.split("_")[3]
    time_df = df["time (min)"]
    mean_df = df.filter(regex="mean")
    std_df = df.filter(regex="std")
    y0 = [500, 500, 0, 0]
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

    params_list = list()
    full_result_dict = dict()

    for condition in conditions:
        if test == "Ratio":
            params_group_list = [condition]
        else:
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

        elif equation == "first_kinetic":
            model = Model(first_kinetic)
            params = model.make_params(k1=dict(value=1e-5, min=10e-9, max=10e1))
            result = model.fit(mean, params, t=time_df, weights=1 / std)
            params_group_list.extend(
                [result.params["k1"].value, result.rsquared, result.params["k1"].stderr]
            )
            if len(labels) != len(params_group_list):
                raise ValueError("the fit output and the equation does not match")

        elif equation == "sec_kinetic":
            model = Model(kin_fit, independent_vars=["t", "y0"])
            params = model.make_params(k1=dict(value=1e-5, min=10e-9, max=10e1))
            result = model.fit(mean, params, t=time_df, y0=y0, weights=1 / std)
            params_group_list.extend(
                [result.params["k1"].value, result.rsquared, result.params["k1"].stderr]
            )
            if len(labels) != len(params_group_list):
                raise ValueError("the fit output and the equation does not match")

        if test != "Ratio":
            result_dict = dict({"{}_{}".format(plate_num, condition): result})
        elif test == "Ratio":
            result_dict = dict({"{}".format(condition): result})

        param_group_dict = dict(zip(labels, params_group_list))
        params_list.append(param_group_dict)
        if not bool(full_result_dict):
            full_result_dict = result_dict
        else:
            full_result_dict.update(result_dict)

    return params_list, full_result_dict


def parameter_determination(test, individual=False):
    one_phase_df_list = list()
    one_phase_result_dict = dict()

    first_kinetic_df_list = list()
    first_kinetic_result_dict = dict()

    sec_kinetic_df_list = list()
    sec_kinetic_result_dict = dict()
    if not individual:
        str_query = "*_{}_*_summerized.pkl".format(test)
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
        str_query = "*_{}_*.pkl".format(test)
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
            first_kinetic_params_list, first_kinetic_results = overall_fit(
                f, test, "first_kinetic", kinetic_list
            )
            sec_kinetic_params_list, sec_kinetic_results = overall_fit(
                f, test, "sec_kinetic", kinetic_list
            )

        else:
            one_phase_params_list, one_phase_results = ind_fit(
                f, test, "one_phase", one_phase_list
            )
            first_kinetic_params_list, first_kinetic_results = ind_fit(
                f, test, "first_kinetic", kinetic_list
            )
            sec_kinetic_params_list, sec_kinetic_results = ind_fit(
                f, test, "sec_kinetic", kinetic_list
            )

        if not bool(one_phase_result_dict):
            one_phase_result_dict = one_phase_results
        else:
            one_phase_result_dict.update(one_phase_results)
        one_phase_df = pd.DataFrame(one_phase_params_list)
        one_phase_df_list.append(one_phase_df)

        if not bool(first_kinetic_result_dict):
            first_kinetic_result_dict = first_kinetic_results
        else:
            first_kinetic_result_dict.update(first_kinetic_results)
        first_kinetic_df = pd.DataFrame(first_kinetic_params_list)
        first_kinetic_df_list.append(first_kinetic_df)

        if not bool(sec_kinetic_result_dict):
            sec_kinetic_result_dict = sec_kinetic_results
        else:
            sec_kinetic_result_dict.update(sec_kinetic_results)
        sec_kinetic_df = pd.DataFrame(sec_kinetic_params_list)
        sec_kinetic_df_list.append(sec_kinetic_df)

    one_phase_params_df = pd.concat(one_phase_df_list)
    first_kinetic_params_df = pd.concat(first_kinetic_df_list)
    sec_kinetic_params_df = pd.concat(sec_kinetic_df_list)

    output = (
        one_phase_params_df,
        first_kinetic_params_df,
        sec_kinetic_params_df,
        one_phase_result_dict,
        first_kinetic_result_dict,
        sec_kinetic_result_dict,
    )

    return output


if __name__ == "__main__":
    individual = True
    if individual:
        os.chdir("./dna_sdr/IO/Output/Pickles/Individual/")
    else:
        os.chdir("./dna_sdr/IO/Output/Pickles/")

    test_type = "Screen"
    output = parameter_determination(test_type, individual=True)

    if individual:
        os.chdir("../../../../..")
    else:
        os.chdir("../../../..")

    os.chdir("./dna_sdr/pickles/")

    if individual:
        if os.path.exists("indiviual_{}_one_phase_param.pkl".format(test_type)):
            one_phase_params_df = pd.read_pickle(
                "indiviual_{}_one_phase_param.pkl".format(test_type)
            )
            if not one_phase_params_df.equals(output[0]):
                output[0].to_pickle(
                    "indiviual_{}_one_phase_param.pkl".format(test_type)
                )
                with open(
                    "indiviual_{}_one_phase_result.pkl".format(test_type), "wb"
                ) as fp:
                    pickle.dump(output[3], fp)
        else:
            output[0].to_pickle("indiviual_{}_one_phase_param.pkl".format(test_type))
            with open(
                "indiviual_{}_one_phase_result.pkl".format(test_type), "wb"
            ) as fp:
                pickle.dump(output[3], fp)

        if os.path.exists("indiviual_{}_first_kinetic_param.pkl".format(test_type)):
            first_kine_params_df = pd.read_pickle(
                "indiviual_{}_first_kinetic_param.pkl".format(test_type)
            )
            if not first_kine_params_df.equals(output[1]):
                output[1].to_pickle(
                    "indiviual_{}_first_kinetic_param.pkl".format(test_type)
                )
                with open(
                    "indiviual_{}_first_kinetic_result.pkl".format(test_type), "wb"
                ) as fp:
                    pickle.dump(output[4], fp)

        else:
            output[1].to_pickle(
                "indiviual_{}_first_kinetic_param.pkl".format(test_type)
            )
            with open(
                "indiviual_{}_first_kinetic_result.pkl".format(test_type), "wb"
            ) as fp:
                pickle.dump(output[4], fp)

        if os.path.exists("indiviual_{}_second_kinetic_param.pkl".format(test_type)):
            sec_kine_params_df = pd.read_pickle(
                "indiviual_{}_second_kinetic_param.pkl".format(test_type)
            )
            if not sec_kine_params_df.equals(output[2]):
                output[2].to_pickle(
                    "indiviual_{}_second_kinetic_param.pkl".format(test_type)
                )
                with open(
                    "indiviual_{}_second_kinetic_result.pkl".format(test_type), "wb"
                ) as fp:
                    pickle.dump(output[5], fp)

        else:
            output[2].to_pickle(
                "indiviual_{}_second_kinetic_param.pkl".format(test_type)
            )
            with open(
                "indiviual_{}_second_kinetic_result.pkl".format(test_type), "wb"
            ) as fp:
                pickle.dump(output[5], fp)

    else:
        if os.path.exists("{}_one_phase_param.pkl".format(test_type)):
            one_phase_params_df = pd.read_pickle(
                "{}_one_phase_param.pkl".format(test_type)
            )
            if not one_phase_params_df.equals(output[0]):
                output[0].to_pickle("{}_one_phase_param.pkl".format(test_type))
                with open("{}_one_phase_result.pkl".format(test_type), "wb") as fp:
                    pickle.dump(output[3], fp)
        else:
            output[0].to_pickle("{}_one_phase_param.pkl".format(test_type))
            with open("{}_one_phase_result.pkl".format(test_type), "wb") as fp:
                pickle.dump(output[3], fp)

        if os.path.exists("{}_first_kinetic_param.pkl".format(test_type)):
            first_kine_params_df = pd.read_pickle(
                "{}_first_kinetic_param.pkl".format(test_type)
            )
            if not first_kine_params_df.equals(output[1]):
                output[2].to_pickle("{}_first_kinetic_param.pkl".format(test_type))
                with open("{}_first_kinetic_result.pkl".format(test_type), "wb") as fp:
                    pickle.dump(output[4], fp)

        else:
            output[2].to_pickle("{}_first_kinetic_param.pkl".format(test_type))
            with open("{}_first_kinetic_result.pkl".format(test_type), "wb") as fp:
                pickle.dump(output[4], fp)

        if os.path.exists("{}_second_kinetic_param.pkl".format(test_type)):
            sec_kine_params_df = pd.read_pickle(
                "{}_second_kinetic_param.pkl".format(test_type)
            )
            if not sec_kine_params_df.equals(output[1]):
                output[2].to_pickle("{}_second_kinetic_param.pkl".format(test_type))
                with open("{}_second_kinetic_result.pkl".format(test_type), "wb") as fp:
                    pickle.dump(output[5], fp)

        else:
            output[2].to_pickle("{}_second_kinetic_param.pkl".format(test_type))
            with open("{}_second_kinetic_result.pkl".format(test_type), "wb") as fp:
                pickle.dump(output[5], fp)
