"""
xxx
"""
import os
import glob
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
from lmfit import Model
from scipy.integrate import solve_ivp
from dna_sdr.experimental.data_processing import time_list_generation
from dna_sdr.visualization.vis import release_grid_plot

plt_params = {
    # "figure.figsize": [6.4, 4.8],
    # "axes.labelsize": 24,
    # "axes.titlesize": 16,
    "axes.linewidth": 2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelpad": 6,
    "axes.formatter.use_mathtext": True,
    "errorbar.capsize": 5,
    "xtick.labelsize": 14,
    "xtick.major.size": 6,
    "xtick.major.width": 2,
    "ytick.labelsize": 14,
    "ytick.major.size": 6,
    "ytick.major.width": 2,
}
pylab.rcParams.update(plt_params)


def one_phase_association(time: float, plateau: float, k: float) -> float:
    """The function calculates the value of a one-phase association reaction over time.

    Parameters
    ----------
    time : float
        The time parameter represents the time at which the association is being calculated. It is a float
    value.
    plateau : float
        The plateau is the maximum value that the function will approach as time goes to infinity.
    k : float
        The parameter "k" represents the rate constant in the one-phase association equation. It determines
    how quickly the association reaction occurs. A higher value of "k" indicates a faster association
    rate.

    Returns
    -------
        the result of the equation `plateau * (1 - np.exp(-k * time))`.

    """
    return plateau * (1 - np.exp(-k * time))


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


def _one_phase_fitting(
    df: pd.DataFrame | pd.Series,
    condition: str,
    labels: list[str],
    individual_fit: bool = False,
):
    """The `_one_phase_fitting` function fits a one-phase association model to data and returns the fitting
    result and a list of parameters.

    Parameters
    ----------
    df : pd.DataFrame | pd.Series
        The `df` parameter is a pandas DataFrame or Series that contains the data for fitting. It is used
    to extract the mean and standard deviation values for fitting the model.
    condition : str
        The `condition` parameter is a string that represents a condition or group in your data. It is used
    to filter and select specific columns from your DataFrame or Series.
    labels : list[str]
        The `labels` parameter is a list of strings that represents the labels for the output parameters.
    These labels will be used to match the fit output with the corresponding equation.
    individual_fit : bool, optional
        A boolean parameter that determines whether to perform an individual fit for each trial or a group
    fit for all trials together. If set to True, the function will fit the model to each trial
    separately. If set to False, the function will fit the model to the mean values of all trials.

    Returns
    -------
        two values: `result` and `params_group_list`.

    """

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
    test_type: str = "Screen",
    individual_fit: bool = False,
    labels: list[str] | None = None,
):
    """The `_kinetic_fitting` function performs kinetic fitting on a given dataset based on the specified
    condition and test type, returning the fitting result and a list of fitting parameters.

    Parameters
    ----------
    df : pd.DataFrame | pd.Series
        The `df` parameter is a pandas DataFrame or Series that contains the data for fitting. It
    represents the experimental data that you want to fit a kinetic model to.
    condition : str
        The `condition` parameter is a string that represents the condition under which the kinetic fitting
    is performed. It is used to filter and process the data accordingly.
    test_type : str, optional
        The `test_type` parameter is a string that specifies the type of test being performed. It can have
    one of three values: "Conc", "Screen", or "Ratio".
    individual_fit : bool, optional
        A boolean flag indicating whether to perform individual fits for each trial or a mean fit for all
    trials combined.
    labels : list[str] | None
        The `labels` parameter is a list of strings that specifies the column labels for the output
    DataFrame or Series. The default value is `None`, in which case the labels will be set to
    `["Condition", "k_rate", "r_sq_kin", "k_error"]`.

    Returns
    -------
        two values: `result` and `params_group_list`.

    """
    if labels is None:
        labels = [
            "Condition",
            "k_rate",
            "r_sq_kin",
            "k_error",
        ]

    if not test_type in {"Conc", "Screen", "Ratio"}:
        raise ValueError()

    if test_type == "Conc":
        factor = int(condition.split("_")[-1]) / 100
        y0 = [500 * factor, 500, 0, 0]
    else:
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
    _type_
        _description_

    """

    df: pd.DataFrame = pd.read_pickle(file)
    name = file.split(".")[0]
    cond = "_".join(name.split("_")[-2:])
    test_type = name.split("_")[2]

    params_list: list[dict[str, str]] = []
    full_result_dict = {}

    for count in range(len(df.columns)):
        trial = df.iloc[:, count] * 500
        if equation == "one_phase":
            result, params_group_list = _one_phase_fitting(trial, cond, labels, True)
        elif equation == "sec_kinetic":
            result, params_group_list = _kinetic_fitting(trial, cond, test_type, True)
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
    """The `overall_fit` function takes in a file, test type, equation type, and labels, and performs
    fitting calculations on the data in the file based on the specified test and equation types,
    returning the fitting parameters and results.

    Parameters
    ----------
    file : str
        The `file` parameter is a string that represents the file name or path of the data file to be read.
    This file should be in pickle format.
    test : str
        The `test` parameter is a string that specifies the type of test being performed. It can have one
    of three values: "Conc", "Screen", or "Ratio".
    equation : str
        The `equation` parameter is a string that specifies the type of equation to use for fitting. It can
    take two possible values: "one_phase" or "sec_kinetic".
    labels : list[str]
        The `labels` parameter is a list of strings that represents the labels for the parameters in the
    fitting equation. These labels are used to create a dictionary that maps each label to its
    corresponding parameter value in the fitting result.

    Returns
    -------
        The function `overall_fit` returns two values: `params_list` and `full_result_dict`.

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
        raise ValueError("Only valid test methods are 'Conc', 'Screen', and 'Ratio'")

    params_list: list[dict[str, str]] = []
    full_result_dict = {}

    for condition in conditions:  # type: ignore - First if statement should ensure it's not unbound
        if equation == "one_phase":
            result, params_group_list = _one_phase_fitting(df, condition, labels)
        elif equation == "sec_kinetic":
            result, params_group_list = _kinetic_fitting(df, condition, test)
        else:
            raise ValueError()

        if test != "Ratio":
            params_group_list.insert(0, plate_num)
            result_dict = dict({f"{plate_num}_{condition}": result})
        elif test == "Ratio":
            result_dict = dict({f"{condition}": result})

        param_group_dict = dict(zip(labels, params_group_list))
        params_list.append(param_group_dict)

        if not full_result_dict:
            full_result_dict = result_dict
        else:
            full_result_dict.update(result_dict)

    return params_list, full_result_dict


# TODO: Fix the dictonary override issue
# * Change the dictonary into a list of dict of dicts with the outer key as the condition and
# * inner key as the trial number
def parameter_determination(
    test: str, individual: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    """_summary_

    Parameters
    ----------
    test : str
        _description_
    individual : bool, optional
        _description_, by default False

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, dict, dict]
        _description_
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
        if not result_dict:
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
    print(sec_kinetic_result_dict)

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
    """The `storage` function saves a DataFrame and a dictionary to pickle files based on the provided
    parameters.

    Parameters
    ----------
    fit_results : tuple[pd.DataFrame, dict]
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
    if individual:
        param_fname = f"individual_{test_type}_{equation}_param.pkl"
        result_fname = f"individual_{test_type}_{equation}_result.pkl"
    else:
        param_fname = f"{test_type}_{equation}_param.pkl"
        result_fname = f"{test_type}_{equation}_result.pkl"

    new_df.to_pickle(param_fname)
    with open(result_fname, "wb") as fp:
        pickle.dump(result_dict, fp)


# TODO: redo the os path finding with os.path.join
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
        os.chdir("./dna_sdr/IO/Output/Individual/")
    else:
        os.chdir("./dna_sdr/IO/Output/Group/")

    output = parameter_determination(test_type, individual)

    os.chdir("../../../..")

    os.chdir("./dna_sdr/pickles/")
    one_phase_results = (output[0], output[2])
    storage(one_phase_results, individual, test_type, "one_phase")

    kinetic_results = (output[1], output[3])
    storage(kinetic_results, individual, test_type, "second_kinetic")


if __name__ == "__main__":
    INDIVIDUAL = True
    TEST_TYPE = "Conc"
    os.chdir("./dna_sdr/IO/Output/Individual/")
    parameter_determination(TEST_TYPE, INDIVIDUAL)

    # main(INDIVIDUAL, TEST_TYPE)

    # df = pd.read_pickle("./dna_sdr/pickles/individual_Conc_second_kinetic_result.pkl")
    # print(df)
    # time = time_list_generation(60)

    # read_lst = ["P1_A1", "P1_A7", "P1_A8", "P1_A9", "P1_A10"]
    # for name in read_lst:
    #     sample: pd.DataFrame = pd.read_pickle(
    #         f"./dna_sdr/IO/Output/Individual/4WJ_HEX_Screen_{name}.pkl"
    #     )
    #     sample = sample * 500
    #     MAX_VALUE = 0
    #     for _, conc in sample.items():
    #         if MAX_VALUE <= conc.max():
    #             MAX_VALUE = conc.max()
    #     y_axis_range = np.arange(0, MAX_VALUE + 50, 100)

    #     fig, axes = plt.subplots(4, 3, figsize=(12, 8), layout="constrained")
    #     ROW = 0
    #     COLUMN = 0

    #     for _, conc in sample.items():
    #         result, _ = _kinetic_fitting(conc, name, individual_fit=True)
    #         release_grid_plot(conc, result, MAX_VALUE, ax=axes[ROW, COLUMN])
    #         if ROW % 3 == 0 and ROW != 0:
    #             COLUMN += 1
    #             ROW = 0
    #         else:
    #             ROW += 1

    #     fig.supylabel("Release quantity (nM)", fontsize=20)
    #     fig.supxlabel("Time (mins)", fontsize=20)
    #     fig.suptitle(f"{name}", fontsize=24)
    #     plt.show()

    # TRIG = "T1"
    # percentage_lst = ["25", "50", "75", "100"]
    # percentage_lst = ["25", "50", "75"]
    # parameters = pd.read_pickle(
    #     "./dna_sdr/pickles/individual_Conc_second_kinetic_param.pkl"
    # )

    # for percentage in percentage_lst:
    #     sample: pd.DataFrame = pd.read_pickle(
    #         f"./dna_sdr/IO/Output/Individual/4WJ_HEX_Conc_{TRIG}_{percentage}.pkl"
    #     )
    #     sample = sample * 500

    #     MAX_VALUE = 0
    #     for _, conc in sample.items():
    #         if MAX_VALUE <= conc.max():
    #             MAX_VALUE = conc.max()
    #     y_axis_range = np.arange(0, MAX_VALUE + 50, 100)

    #     fig, axes = plt.subplots(4, 3, figsize=(12, 8), layout="constrained")
    #     ROW = 0
    #     COLUMN = 0

    #     for _, conc in sample.items():
    #         print()
    #         release_grid_plot(conc, result, MAX_VALUE, ax=axes[ROW, COLUMN])
    #         if ROW % 3 == 0 and ROW != 0:
    #             COLUMN += 1
    #             ROW = 0
    #         else:
    #             ROW += 1

    #     fig.supylabel("Release quantity (nM)", fontsize=20)
    #     fig.supxlabel("Time (mins)", fontsize=20)
    #     fig.suptitle(f"{TRIG}_{percentage}", fontsize=24)
    #     plt.show()
