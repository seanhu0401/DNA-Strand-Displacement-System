import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from lmfit import minimize, Parameters


def simple(t, X, const):
    # simplified system of ODE for SDR system
    I = X[0]
    QT = X[1]
    IT = X[2]
    Q = X[3]

    try:
        k1 = const["k1"].value

    except KeyError:
        k1 = const

    # the model equations
    dIdt = -k1 * I * QT
    dQTdt = -k1 * I * QT
    dITdt = k1 * I * QT
    dQdt = k1 * I * QT
    return [dIdt, dQTdt, dITdt, dQdt]


def fn_solve(t, x0, paras):
    x = solve_ivp(
        simple, (0, max(t) + 10), x0, t_eval=t, args=(paras,), rtol=1e-9, method="LSODA"
    )
    return x.y


def residual(paras, t, data):
    """
    compute the residual between actual data and fitted data
    """
    C0 = (
        paras["I"].value,
        paras["QT"].value,
        paras["IT"].value,
        paras["Q"].value,
    )
    model = fn_solve(t, C0, paras)

    IQ_model = model[2]
    return (IQ_model - data).ravel()


def IVP_fitting(file):
    df = pd.read_pickle(file)
    plate_num = file.split("_")[3]
    time = df["time (min)"]
    mean_df = df.filter(regex="mean")
    conditions = [i.split("_")[0] for i in mean_df.columns]

    master_list = list()
    label_list = [
        "Plate Number",
        "Trigger type",
        "rate",
        "r_sq",
        "rate_error",
    ]

    for condition in conditions:
        y0 = np.zeros(4)
        y0[0:2] = 500

        params = Parameters()
        params.add("I", value=y0[0], vary=False)
        params.add("QT", value=y0[1], vary=False)
        params.add("IT", value=y0[2], vary=False)
        params.add("Q", value=y0[3], vary=False)
        params.add("k1", value=1e-6, min=10e-9, max=10e1)

        group_list = [plate_num, condition]
        mean = mean_df.filter(regex=condition) * 500
        mean = mean.squeeze()
        result = minimize(residual, params, args=(time, mean), method="leastsq")
        r2 = 1 - result.residual.var() / np.var(mean)

        group_list.extend([result.params["k1"].value, r2, result.params["k1"].stderr])
        param_group_dict = dict(zip(label_list, group_list))
        master_list.append(param_group_dict)

    return master_list


if __name__ == "__main__":
    os.chdir("./dna_sdr/IO/Output/Pickles")
    df_list = list()
    for f in glob.glob("*" + "Screen" + "*" + "summerized.pkl"):
        rate_list = IVP_fitting(f)
        df = pd.DataFrame(rate_list)
        df_list.append(df)

    parameter_df = pd.concat(df_list)
    parameter_df = parameter_df.mask(parameter_df["r_sq"] <= 0.80).dropna()

    trig_params_df = parameter_df.loc[parameter_df["Trigger type"] != "T1"]

    # print(trig_params_df)
    trig_params_df.to_pickle("trig_kinetic.pkl")
    trig_params_df.to_csv("trig_kinetic.csv")
