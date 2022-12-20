import pandas as pd
import numpy as np
import warnings
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")


def general_cf_process(eq, x_df, y_df):
    r_square = 0
    popt = None
    try:
        popt, pcov = curve_fit(eq, x_df, y_df)
        residual = y_df - eq(np.array(x_df), *popt)
        ss_res = np.sum(residual**2)
        ss_tot = np.sum((y_df - np.mean(y_df)) ** 2)
        r_square = 1 - (ss_res / ss_tot)

    except RuntimeError:
        pass

    return popt, r_square
