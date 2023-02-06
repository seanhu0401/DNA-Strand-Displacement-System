import pandas as pd
import numpy as np
import warnings
from scipy.optimize import curve_fit
from dna_sdr.curve_fitting.exp_models import one_phase_association

warnings.filterwarnings("ignore")


def general_cf_process(eq, x_df, y_df):
    r_square = 0
    popt = None
    try:
        popt, _ = curve_fit(eq, x_df, y_df)
        residual = y_df - eq(np.array(x_df), *popt)
        ss_res = np.sum(residual**2)
        ss_tot = np.sum((y_df - np.mean(y_df)) ** 2)
        r_square = 1 - (ss_res / ss_tot)

    except RuntimeError:
        pass

    return popt, r_square


def Individual_Curve_Fit(norm_release, time_df):
    master_list = list()
    column_name = ["sample name", "plateau", "rate (min^-1)", "r-square"]

    for columns in norm_release:
        try:
            popt, r_square = general_cf_process(
                one_phase_association, time_df, norm_release[columns]
            )
            if r_square > 0.8:
                popt_list = popt.tolist()
                popt_list.insert(0, columns)
                popt_list.append(r_square)
                master_list.append(popt_list)
        except RuntimeError:
            pass

    master_df = pd.DataFrame(master_list, columns=column_name)
    sorted_master_df = master_df.sort_values(by=["sample name"])

    return sorted_master_df


def Average_Curve_Fit(exp_release_df):
    time_df = exp_release_df["time (min)"]
    norm_release_df = exp_release_df.filter(regex="mean")

    variable_list = list()
    fitted_value_list = list()
    fitted_value_col_name = list()
    column_name = ["Type", "Plateau", "k", "r-square"]

    for column in norm_release_df:
        avg_release = norm_release_df[column]
        split_col = column.split("_")
        col_name = f"{split_col[0]} ({split_col[1]})"
        try:
            popt, r_square = general_cf_process(
                one_phase_association, time_df, avg_release
            )
            fitted_value = one_phase_association(np.array(time_df), *popt)
            popt_list = popt.tolist()
            popt_list.insert(0, col_name)
            popt_list.append(r_square)
            variable_list.append(popt_list)
            fitted_value_list.append(fitted_value)
            fitted_value_col_name.append(col_name)
        except RuntimeError:
            continue

    variable_df = pd.DataFrame(variable_list, columns=column_name)
    fitted_value_df = pd.DataFrame(
        fitted_value_list, columns=time_df, index=fitted_value_col_name
    ).T.reset_index()
    fitted_value_df.set_index("time (min)", inplace=True)

    return variable_df, fitted_value_df
