import glob
import os
import pandas as pd
import dna_sdr.fitting.curve as cf
import dna_sdr.fitting.kinetic as kinetic


def param_df_gen():
    curve_df_list = list()
    curve_list = [
        "Plate Number",
        "Trigger type",
        "plateau",
        "rate",
        "r_sq_curve",
        "plateau_error",
        "rate_error",
    ]

    kinetic_df_list = list()
    kinetic_list = [
        "Plate Number",
        "Trigger type",
        "k_rate",
        "r_sq_kin",
        "k_error",
    ]

    for f in glob.glob("*" + "Screen" + "*" + "summerized.pkl"):
        curve_parameter_list = cf.cf(f, curve_list)
        curve_df = pd.DataFrame(curve_parameter_list)
        curve_df_list.append(curve_df)

        kinetic_parameter_list = kinetic.IVP_fitting(f, kinetic_list)
        kinetic_df = pd.DataFrame(kinetic_parameter_list)
        kinetic_df_list.append(kinetic_df)

    curve_parameter_df = pd.concat(curve_df_list)
    kinetic_parameter_df = pd.concat(kinetic_df_list)

    return curve_parameter_df, kinetic_parameter_df


if __name__ == "__main__":
    os.chdir("./dna_sdr/IO/Output/Pickles")
    curve_param_df, kinetic_param_df = param_df_gen()

    trig_curve_param_df = curve_param_df.loc[curve_param_df["Trigger type"] != "T1"]
    trig_kin_param_df = kinetic_param_df.loc[kinetic_param_df["Trigger type"] != "T1"]

    os.chdir("../../../..")
    trig_curve_param_df.to_pickle("./dna_sdr/pickles/trig_curve_param.pkl")
    trig_kin_param_df.to_pickle("./dna_sdr/pickles/trig_kin_param.pkl")

    trig_curve_param_df["Plate loc"] = (
        trig_curve_param_df["Plate Number"] + "_" + trig_curve_param_df["Trigger type"]
    )
    trig_curve_param_df.drop(["Plate Number", "Trigger type"], axis=1, inplace=True)

    trig_kin_param_df["Plate loc"] = (
        trig_kin_param_df["Plate Number"] + "_" + trig_kin_param_df["Trigger type"]
    )
    trig_kin_param_df.drop(["Plate Number", "Trigger type"], axis=1, inplace=True)

    param_df = pd.merge(trig_curve_param_df, trig_kin_param_df, on="Plate loc")
    param_df.to_pickle("./dna_sdr/pickles/trig_param.pkl")
