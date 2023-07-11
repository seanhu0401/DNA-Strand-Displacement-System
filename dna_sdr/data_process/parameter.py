import os
import pandas as pd

if __name__ == "__main__":
    os.chdir("./dna_sdr/pickles")
    curve_param_df = pd.read_pickle("screen_one_phase_param.pkl")
    kinetic_param_df = pd.read_pickle("screen_kinetic_param.pkl")

    trig_curve_param_df = curve_param_df.loc[curve_param_df["Trigger type"] != "T1"]
    trig_kin_param_df = kinetic_param_df.loc[kinetic_param_df["Trigger type"] != "T1"]

    trig_curve_param_df.to_pickle("trig_one_phase_param.pkl")
    trig_kin_param_df.to_pickle("trig_kin_param.pkl")

    trig_curve_param_df["plate_loc"] = (
        trig_curve_param_df["Plate Number"] + "_" + trig_curve_param_df["Trigger type"]
    )
    trig_curve_param_df.drop(["Plate Number", "Trigger type"], axis=1, inplace=True)

    trig_kin_param_df["plate_loc"] = (
        trig_kin_param_df["Plate Number"] + "_" + trig_kin_param_df["Trigger type"]
    )
    trig_kin_param_df.drop(["Plate Number", "Trigger type"], axis=1, inplace=True)

    param_df = pd.merge(trig_curve_param_df, trig_kin_param_df, on="plate_loc")
    param_df.to_pickle("trig_param.pkl")
