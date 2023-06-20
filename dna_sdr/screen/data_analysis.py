import glob
import os
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit


def one_phase_association(time, plateau, k):
    return plateau * (1 - np.exp(-k * time))


def general_cf_process(eq, x_df, y_df, stdev):
    r_square = 0
    popt = [0, 0]
    perr = [0, 0]
    try:
        popt, pcov = curve_fit(
            eq, x_df, y_df, sigma=stdev, absolute_sigma=True, p0=[500, 0.01]
        )
        residual = y_df - eq(np.array(x_df), *popt)
        ss_res = np.sum(residual**2)
        ss_tot = np.sum((y_df - np.mean(y_df)) ** 2)
        r_square = 1 - (ss_res / ss_tot)
        perr = np.sqrt(np.diag(pcov))

    except RuntimeError:
        pass

    return (popt, r_square, perr)


def CF(file):
    if not os.path.exists(file):
        os.chdir("./dna_sdr/IO/Output/Pickles")
    print(file)

    df = pd.read_pickle(file)
    plate_num = file.split("_")[3]
    time = df["time (min)"]
    mean_df = df.filter(regex="mean")
    std_df = df.filter(regex="std")
    conditions = [i.split("_")[0] for i in mean_df.columns]

    master_list = list()
    label_list = [
        "Plate Number",
        "Trigger type",
        "plateau",
        "rate",
        "r_sq",
        "plateau_error",
        "rate_error",
    ]

    for condition in conditions:
        group_list = [plate_num, condition]
        mean = mean_df.filter(regex=condition) * 500
        std = std_df.filter(regex=condition) * 500
        mean = mean.squeeze()
        std = std.squeeze()
        results = general_cf_process(one_phase_association, time, mean, std)

        group_list.extend([*results[0], results[1], *results[2]])

        param_group_dict = dict(zip(label_list, group_list))
        master_list.append(param_group_dict)

    return master_list


def trig_param_df_gen():
    df_list = list()
    for f in glob.glob("*" + "Screen" + "*" + "summerized.pkl"):
        parameter_list = CF(f)
        df = pd.DataFrame(parameter_list)
        df_list.append(df)
        parameter_df = pd.concat(df_list)
        parameter_df = parameter_df.mask(parameter_df["r_sq"] <= 0.90).dropna()

    trig_param_df = parameter_df.loc[parameter_df["Trigger type"] != "T1"]

    return trig_param_df


if __name__ == "__main__":
    os.chdir("./dna_sdr/IO/Output/Pickles")
    trig_df = trig_param_df_gen()
    print(trig_df)

    # f = "4WJ_HEX_Screen_P0_A01_A04_summerized.pkl"
    # parameter_list = CF(f)
    # df = pd.DataFrame(parameter_list)
    # print(df)

    # os.chdir("../../../..")
    # trig_df.to_pickle("./dna_sdr/pickles/trig.pkl")
    # trig_df.to_csv("./dna_sdr/trig-df.csv")
