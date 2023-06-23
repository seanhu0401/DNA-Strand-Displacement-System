import pandas as pd
import numpy as np
import warnings
import os
from scipy.optimize import curve_fit


warnings.filterwarnings("ignore")

"""
The equation(s) below was used to fit the release curve from the qPCR experiments.
"""


def one_phase_association(time, plateau, k):
    # pseudo-first order association kinetics
    return plateau * (1 - np.exp(-k * time))


def general_cf_process(eq, x_df, y_df):
    r_square = 0
    popt = None
    perr = None
    try:
        popt, pcov = curve_fit(eq, x_df, y_df)
        residual = y_df - eq(np.array(x_df), *popt)
        ss_res = np.sum(residual**2)
        ss_tot = np.sum((y_df - np.mean(y_df)) ** 2)
        r_square = 1 - (ss_res / ss_tot)
        perr = np.sqrt(np.diag(pcov))

    except RuntimeError:
        pass

    return (popt, r_square, perr)


def cf(file, labels: list, eq=one_phase_association):
    master_list = list()
    if not os.path.exists(file):
        os.chdir("./dna_sdr/IO/Output/Pickles")
    print(file)

    df = pd.read_pickle(file)
    plate_num = file.split("_")[3]
    time = df["time (min)"]
    mean_df = df.filter(regex="mean")
    std_df = df.filter(regex="std")
    conditions = [i.split("_")[0] for i in mean_df.columns]

    for condition in conditions:
        group_list = [plate_num, condition]
        mean = mean_df.filter(regex=condition) * 500
        std = std_df.filter(regex=condition) * 500
        mean = mean.squeeze()
        std = std.squeeze()
        results = general_cf_process(eq, time, mean, std)

        group_list.extend([*results[0], results[1], *results[2]])

        if len(group_list) != len(labels):
            print(labels)
            print("length of parameter list is {}".format(len(group_list)))
            print("length of label list is {}".format(len(labels)))
            raise ValueError(
                "Length of the list of labels does not match the length of the parameters"
            )
        else:
            param_group_dict = dict(zip(labels, group_list))

        master_list.append(param_group_dict)
    return master_list
