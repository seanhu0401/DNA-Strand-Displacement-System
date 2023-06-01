import glob
import os
import pandas as pd
import numpy as np


def name_conversion(trig_df, info_df):
    plate_loc_lst = list()
    name_lst = list()

    for index in trig_df.index:
        plate_loc = "_".join(index)
        info = info_df[info_df["plate_loc"] == plate_loc]
        name = info["name"]
        if not name.empty:
            plate_loc_lst.append(plate_loc)
            name_lst.append(name.iloc[0])

    convert_dict = dict(zip(plate_loc_lst, name_lst))
    return convert_dict


if __name__ == "__main__":
    t1_pickle = "./dna_sdr/pickles/T1_pt.pkl"
    trig_pickle = "./dna_sdr/pickles/trig_pt.pkl"
    conc_pickle = "./dna_sdr/pickles/concentration.pkl"
    trig_info_pickle = "./dna_sdr/pickles/trig_info.pkl"

    t1_df = pd.read_pickle(t1_pickle)
    trig_df = pd.read_pickle(trig_pickle)
    trig_info_df = pd.read_pickle(trig_info_pickle)
    conc_df = pd.read_pickle(conc_pickle)

    plateau_series = trig_df["mean"]["plateau"] * 500
    name_convert_dict = name_conversion(trig_df, trig_info_df)
    plate_loc_lst = list()

    for index in plateau_series.index:
        plate_loc = "_".join(index)
        plate_loc_lst.append(plate_loc)

    plateau_df = plateau_series.reset_index()
    plateau_df["plate_loc"] = plate_loc_lst
    plateau_df.drop(["Plate Number", "Trigger type"], axis=1, inplace=True)
    trig_info_plateau_df = pd.merge(trig_info_df, plateau_df, on="plate_loc")

    print(trig_info_plateau_df)
    print(conc_df)
