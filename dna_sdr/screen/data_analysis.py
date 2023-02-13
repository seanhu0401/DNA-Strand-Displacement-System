import glob
import os
import pandas as pd
import numpy as np
import dna_sdr.data_process.list_generation as lst_gen
import dna_sdr.data_process.group as grouping
import dna_sdr.curve_fitting.curve_fitting as cf

# TODO: Change the file directory for better os paths


def Individual_CF_Parameter(file):
    if not os.path.exists(file):
        os.chdir("./dna_sdr/IO/Output/Pickles")
    df = pd.read_pickle(file)
    plate_number, group_dict = grouping.screen_grouping(file)
    time_df = lst_gen.time_list_generation(len(df))
    master_list = list()
    label_list = [
        "Plate Number",
        "Trigger type",
        "plateau",
        "rate",
        "r_sq",
    ]

    for group in group_dict:
        group_list = [plate_number, group]

        plateau_list, rate_const_list, r_square_list = (list() for i in range(3))
        cond_group = df[group_dict[group]]
        for trial in range(cond_group.shape[1]):
            ind_trial = cond_group.iloc[:, trial]
            parameter, r_square = cf.general_cf_process(
                cf.one_phase_association, time_df, ind_trial
            )
            if parameter is not None:
                plateau_list.append(parameter[0])
                rate_const_list.append(parameter[1])
                r_square_list.append(r_square)
        group_list.append(plateau_list)
        group_list.append(rate_const_list)
        group_list.append(r_square_list)
        param_group_dict = dict(zip(label_list, group_list))
        master_list.append(param_group_dict)

    return master_list


def parameter_df_gen(params_list):
    df_list = list()
    for condition in range(len(params_list)):
        ind_cond_df = pd.DataFrame.from_dict(params_list[condition])
        df_list.append(ind_cond_df)
        parameter_df = pd.concat(df_list)
        parameter_df = parameter_df.mask(parameter_df["r_sq"] <= 0.90).dropna()

    return parameter_df


def pt_gen(ext="pkl"):
    params_df_list = list()
    os.chdir("./dna_sdr/IO/Output/Pickles")
    for f in glob.glob("*" + "Screen" + "*" + "normalized.{}".format(ext)):
        print(f)
        individual_parameter_list = Individual_CF_Parameter(f)
        params_df = parameter_df_gen(individual_parameter_list)
        params_df_list.append(params_df)

    parameters_df = pd.concat(params_df_list)

    T1_df = parameters_df.loc[parameters_df["Trigger type"] == "T1"]

    T1_pt = pd.pivot_table(
        data=T1_df,
        index=["Plate Number", "Trigger type"],
        values=["plateau", "rate"],
        aggfunc=[np.mean, np.std, "count"],
    )

    pt_output = pd.pivot_table(
        data=parameters_df.loc[parameters_df["Trigger type"] != "T1"],
        index=["Plate Number", "Trigger type"],
        values=["plateau", "rate"],
        aggfunc=[np.mean, np.std, "count"],
    )

    return T1_pt, pt_output


if __name__ == "__main__":
    t1_out, trig_out = pt_gen()
    t1_pickle = "./dna_sdr/pickles/T1_pt.pkl"
    trig_pickle = "./dna_sdr/pickles/trig_pt.pkl"

    os.chdir("../../../..")
    if os.path.exists(t1_pickle) and os.path.exists(trig_pickle):
        t1_df = pd.read_pickle(t1_pickle)
        trig_df = pd.read_pickle(trig_pickle)
        if t1_out.equals(t1_df):
            print("Same T1 dataframe - no new T1 file")
        else:
            t1_out.to_pickle(t1_pickle)

        if trig_out.equals(trig_df):
            print("Same trigger dataframe - no new trigger file")
        else:
            trig_out.to_pickle(trig_pickle)
    else:
        t1_out.to_pickle(t1_pickle)
        trig_out.to_pickle(trig_pickle)
