import glob
import os
import pandas as pd
import dna_sdr.fitting.curve as cf


def param_df_gen():
    df_list = list()
    label_list = [
        "Plate Number",
        "Trigger type",
        "plateau",
        "rate",
        "r_sq",
        "plateau_error",
        "rate_error",
    ]

    for f in glob.glob("*" + "Screen" + "*" + "summerized.pkl"):
        parameter_list = cf.cf(f, label_list)
        df = pd.DataFrame(parameter_list)
        df_list.append(df)

    parameter_df = pd.concat(df_list)
    parameter_df = parameter_df.mask(parameter_df["r_sq"] <= 0.90).dropna()

    return parameter_df


if __name__ == "__main__":
    os.chdir("./dna_sdr/IO/Output/Pickles")
    param_df = param_df_gen()
    trig_param_df = param_df.loc[param_df["Trigger type"] != "T1"]

    os.chdir("../../../..")
    trig_param_df.to_pickle("./dna_sdr/pickles/trig.pkl")
    # trig_df.to_csv("./dna_sdr/trig-df.csv")
