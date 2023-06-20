import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from lmfit import minimize, Parameters, report_fit


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
    # print(sum((IQ_model - data).ravel()))
    return (IQ_model - data).ravel()


if __name__ == "__main__":
    os.chdir("./dna_sdr/IO/Output/Pickles")
    df_list = list()
    for f in glob.glob("*" + "Screen" + "*" + "summerized.pkl"):
        df = pd.read_pickle(f)
        plate_num = f.split("_")[3]
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

            simp_params = Parameters()
            simp_params.add("I", value=y0[0], vary=False)
            simp_params.add("QT", value=y0[1], vary=False)
            simp_params.add("IT", value=y0[2], vary=False)
            simp_params.add("Q", value=y0[3], vary=False)
            simp_params.add("k1", value=1e-6, min=10e-9, max=10e-4)

            group_list = [plate_num, condition]
            mean = mean_df.filter(regex=condition) * 500
            mean = mean.squeeze()
            result = minimize(
                residual, simp_params, args=(time, mean), method="leastsq"
            )
            r2 = 1 - result.residual.var() / np.var(mean)

            (
                rate_const_list,
                r_square_list,
                rate_const_error_list,
            ) = (list() for i in range(3))

            rate_const_list.append(result.params["k1"].value)
            r_square_list.append(r2)
            rate_const_error_list.append(result.params["k1"].stderr)

            group_list.append(rate_const_list)
            group_list.append(r_square_list)
            group_list.append(rate_const_error_list)

            param_group_dict = dict(zip(label_list, group_list))
            master_list.append(param_group_dict)

            # print(result.params["k1"].value)
            # print(result.params["k1"].stderr)
            # print("r^2 value is {}".format(r2))
        for condition in range(len(master_list)):
            df = pd.DataFrame.from_dict(master_list[condition])
            df_list.append(df)
        parameter_df = pd.concat(df_list)
        parameter_df = parameter_df.mask(parameter_df["r_sq"] <= 0.90).dropna()
    print(parameter_df)
    # parameter_df.to_csv("kinetic.csv")
