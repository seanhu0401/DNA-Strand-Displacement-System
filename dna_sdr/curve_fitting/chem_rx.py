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


def simple_complex(t, X, const):
    """
    more complex system of ODE for SDR system that includes the formation of intermediate (IQT).
    """

    I = X[0]
    QT = X[1]
    IT = X[2]
    Q = X[3]
    IQT = X[4]

    try:
        k1 = const["k1"].value
        # k1r = const["k1r"].value
        k2 = const["k2"].value

    except KeyError:
        # k1, k1r, k2 = const
        k1, k2 = const

    dIdt = -k1 * QT * I
    dQTdt = -k1 * QT * I
    dIQTdt = k1 * I * QT - k2 * IQT
    dITdt = k2 * IQT
    dQdt = k2 * IQT

    return [dIdt, dQTdt, dITdt, dQdt, dIQTdt]


def complex(t, X, const):
    """
    more complex system of ODE for SDR system that includes the formation of intermediate (IQT).
    """

    I = X[0]
    QT = X[1]
    IT = X[2]
    Q = X[3]
    IQT = X[4]

    try:
        k1 = const["k1"].value
        k1r = const["k1r"].value
        k2 = const["k2"].value

    except KeyError:
        k1, k1r, k2 = const
        k1, k2 = const

    dIdt = -k1 * QT * I + k1r * IQT
    dQTdt = -k1 * QT * I + k1r * IQT
    dIQTdt = k1 * I * QT - k2 * IQT
    dITdt = k2 * IQT
    dQdt = k2 * IQT

    return [dIdt, dQTdt, dITdt, dQdt, dIQTdt]


def fn_solve(t, x0, paras, fn):
    x = solve_ivp(
        fn, (0, max(t) + 10), x0, t_eval=t, args=(paras,), rtol=1e-9, method="LSODA"
    )
    return x.y


def residual(paras, t, data, mode, fn):
    """
    compute the residual between actual data and fitted data
    """
    if mode == "simple":
        C0 = (
            paras["I"].value,
            paras["QT"].value,
            paras["IT"].value,
            paras["Q"].value,
        )
    else:
        C0 = (
            paras["I"].value,
            paras["QT"].value,
            paras["IQT"].value,
            paras["IT"].value,
            paras["Q"].value,
        )

    model = fn_solve(t, C0, paras, fn)

    IQ_model = model[2]
    print(sum((IQ_model - data).ravel()))
    return (IQ_model - data).ravel()


# TODO: Update the graphs to show the model for trig+output in one subfigure + model output with exp output in another subfigure
if __name__ == "__main__":
    pickle = "dna_sdr/IO/Output/Pickles/4WJ_HEX_Screen_P0_A01_A04_summerized.pkl"
    df = pd.read_pickle(pickle)

    # measured data
    t_measured = df["time (min)"]
    release_measured = df["A3_mean"] * 500
    y0 = np.zeros(5)
    y0[0:2] = 500
    k1 = 1e-6
    k1r = 1e-5
    k2 = 1e-0
    k = [k1, k1r, k2]

    k1_bound = (10e-9, 10e-4)
    k1r_bound = (10e-9, 10e-4)
    k2_bound = (10e-9, 10e1)
    k_bound = [k1_bound, k1r_bound, k2_bound]

    # mode = "simple"
    mode = "simple-complex"
    # mode = "complex"

    if mode == "simple":
        # initial conditions
        simp_params = Parameters()
        simp_params.add("I", value=y0[0], vary=False)
        simp_params.add("QT", value=y0[1], vary=False)
        simp_params.add("IT", value=y0[2], vary=False)
        simp_params.add("Q", value=y0[3], vary=False)
        simp_params.add("k1", value=k[0], min=k_bound[0][0], max=k_bound[0][1])

        result = minimize(
            residual,
            simp_params,
            args=(t_measured, release_measured, mode, simple),
            method="leastsq",
        )  # leastsq nelder

        # check results of the fit
        data_fitted = fn_solve(t_measured, y0[:4], result.params, simple)

        # plot fitted data
        plt.scatter(
            t_measured,
            release_measured,
            marker="o",
            color="b",
            label="exp",
            s=75,
        )
        plt.plot(
            t_measured,
            data_fitted[0],
            "-",
            linewidth=2,
            color="black",
            label="I",
        )
        plt.plot(
            t_measured,
            data_fitted[2],
            "-",
            linewidth=2,
            color="red",
            label="QT",
        )
        plt.legend()
        plt.xlim([0, max(t_measured)])
        plt.ylim([0, 1.1 * 500])

        # display fitted statistics
        report_fit(result)
        r2 = 1 - result.residual.var() / np.var(release_measured)
        print("r^2 value is {}".format(r2))

        plt.show()

    elif mode == "simple-complex":
        # initial conditions
        simple_complex_params = Parameters()
        simple_complex_params.add("I", value=y0[0], vary=False)
        simple_complex_params.add("QT", value=y0[1], vary=False)
        simple_complex_params.add("IT", value=y0[2], vary=False)
        simple_complex_params.add("Q", value=y0[3], vary=False)
        simple_complex_params.add("IQT", value=y0[4], vary=False)

        simple_complex_params.add(
            "k1", value=k[0], min=k_bound[0][0], max=k_bound[0][1]
        )
        simple_complex_params.add(
            "k2", value=k[2], min=k_bound[2][0], max=k_bound[2][1]
        )

        result = minimize(
            residual,
            simple_complex_params,
            args=(t_measured, release_measured, mode, simple_complex),
            method="leastsq",
        )  # leastsq nelder

        data_fitted = fn_solve(t_measured, y0, result.params, simple_complex)

        # plot fitted data
        plt.scatter(
            t_measured,
            release_measured,
            marker="o",
            color="b",
            label="measured data",
            s=75,
        )
        plt.plot(
            t_measured,
            data_fitted[0],
            "-",
            linewidth=2,
            color="black",
            label="I",
        )
        plt.plot(
            t_measured,
            data_fitted[2],
            "-",
            linewidth=2,
            color="red",
            label="QT",
        )
        plt.plot(
            t_measured,
            data_fitted[4],
            "-",
            linewidth=2,
            color="green",
            label="IQT",
        )
        plt.legend()
        plt.xlim([0, max(t_measured)])
        plt.ylim([0, 1.1 * 500])

        # display fitted statistics
        report_fit(result)
        r2 = 1 - result.residual.var() / np.var(release_measured)
        print("r^2 value is {}".format(r2))

        plt.show()

    elif mode == "complex":
        # initial conditions
        complex_params = Parameters()
        complex_params.add("I", value=y0[0], vary=False)
        complex_params.add("QT", value=y0[1], vary=False)
        complex_params.add("IT", value=y0[2], vary=False)
        complex_params.add("Q", value=y0[3], vary=False)
        complex_params.add("IQT", value=y0[4], vary=False)

        complex_params.add("k1", value=k[0], min=k_bound[0][0], max=k_bound[0][1])
        complex_params.add("k1r", value=k[1], min=k_bound[1][0], max=k_bound[1][1])
        complex_params.add("k2", value=k[2], min=k_bound[2][0], max=k_bound[2][1])

        result = minimize(
            residual,
            complex_params,
            args=(t_measured, release_measured, mode, complex),
            method="leastsq",
        )  # leastsq nelder
        data_fitted = fn_solve(t_measured, y0, result.params, complex)

        # plot fitted data
        plt.scatter(
            t_measured,
            release_measured,
            marker="o",
            color="b",
            label="exp",
            s=75,
        )
        plt.plot(
            t_measured,
            data_fitted[0],
            "-",
            linewidth=2,
            color="black",
            label="QT",
        )
        plt.plot(
            t_measured,
            data_fitted[2],
            "-",
            linewidth=2,
            color="red",
            label="QT",
        )
        plt.plot(
            t_measured,
            data_fitted[4],
            "-",
            linewidth=2,
            color="green",
            label="IQT",
        )
        plt.legend()
        plt.xlim([0, max(t_measured)])
        plt.ylim([0, 1.1 * 500])

        # display fitted statistics
        report_fit(result)
        r2 = 1 - result.residual.var() / np.var(release_measured)
        print("r^2 value is {}".format(r2))

        plt.show()
