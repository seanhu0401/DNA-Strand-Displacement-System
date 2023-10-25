"""
This module is used to generate the various DataFrame for data analysis

The main purpose of this module is to define the functions that is used to 
process the data obtained from the 'Q' qPCR machine and generate various DataFrame from it.
It combines, normalizes, averages, summerizes the data collected. 

"""

import os
import glob
import regex as re
import pandas as pd


def trig_list_gen(fname_list: list[str]) -> list[str]:
    """
    Generate a list of triggers that are presented within the specific test file.

    Parameters
    ----------
    fname_list : list[str]
        A list that contained all the file name for a specific conditon

    Returns
    -------
    trig_list : list[str]
        A list that contain all the trigger name presented within the specific test file.

    """
    trig_list = list(str())
    re_query = r"[a-zA-Z]\d\d"
    for fname in fname_list:
        split_name = fname.split(".")[0].split("_")
        match = re.findall(re_query, fname)

        if split_name[2] != "Ratio":
            plate_number = split_name[3]
            trig_type = plate_number + "_" + "_".join(match)

            if trig_type not in trig_list:
                trig_list.append(trig_type)

        elif split_name[2] == "Ratio":
            plate_number_1 = split_name[3]
            plate_number_2 = split_name[5]
            trig_type = (
                plate_number_1 + "_" + match[0] + "_" + plate_number_2 + "_" + match[1]
            )

            if trig_type not in trig_list:
                trig_list.append(trig_type)

    return trig_list


def file_path_generation(
    test_type: str, trig_type: str
) -> tuple[str, str, str, str, str]:
    """
    Generate the file paths that are used to save the various dataframe into pickle files.

    Parameters
    ----------
    test_type : str
        The type of test that is being analyized

    trig_type : str
        The triggers that are selected in the file

    Returns
    -------
    dir_name : str
        Name of the save directory

    cdata_path : str
        The sava path for the combined data

    normalized_cdata_path : str
        The save path for the normalized combined data

    summerized_data_path : str
        The save path for the summerized data

    dir_path : str
        The save path for the raw data after processing.

    """
    pickle_extension = "pkl"
    pickle_save_path = "../Output/Pickles"
    parent_dir = "../Processed/"

    if test_type in {"Ratio", "Screen", "Conc"}:
        dir_name = f"4WJ_HEX_{test_type}_{trig_type}"
    else:
        raise ValueError(f"unknown test type - {test_type}")

    combined_data_after_baseline_subtraction_pickle = f"{dir_name}.{pickle_extension}"
    norm_combined_data_pickle = f"{dir_name}_normalized.{pickle_extension}"
    summerized_data_pickle = f"{dir_name}_summerized.{pickle_extension}"
    cdata_path = os.path.join(
        pickle_save_path, combined_data_after_baseline_subtraction_pickle
    )
    normalized_cdata_path = os.path.join(pickle_save_path, norm_combined_data_pickle)
    summarized_data_path = os.path.join(pickle_save_path, summerized_data_pickle)
    dir_path = os.path.join(parent_dir, dir_name)

    return (
        dir_name,
        cdata_path,
        normalized_cdata_path,
        summarized_data_path,
        dir_path,
    )


def time_list_generation(num_cycle: int, init_time_sec: int = 231) -> list[float]:
    """
    Generate a list of time point for the study

    Parameters
    ----------
    num_cycle : int
        The number of cycles that was carried out in the study

    init_time_sec : int = 231
        The starting time of the first cycle. Default to be 231 second due
        to the initial lag time of the machine

    Returns
    -------
    t_list_mins : list[float]
        The list of time points at each measurement time (end of each cycle) in mintues.

    """
    cycle_list = list(range(1, num_cycle))
    t_list_seconds = [init_time_sec]

    for cycle in cycle_list:
        t = t_list_seconds[-1] + 20 + 5 * cycle
        t_list_seconds.append(t)

    t_list_mins = [round(t / 60, 2) for t in t_list_seconds]

    return t_list_mins


def group_query(fname: str, option: str | None = None) -> tuple[int, list[str]]:
    """
    Search the testing groups in the file based on the file name. Calculate how
    many groups are presented in the test and what are the groups.

    Parameters
    ----------
    fname : str
        The name of the test file.

    option: str or None
        An optional augment used for ratio study to determine the if the group is
        over or under 100%. The default value it 'None'. The valid option input
        strings are 'over' or 'under'.

    Returns
    -------
    group_quant : int
        The number of groups presented in the test file.

    group_lst: list[str]
        The groups presented in the test file.

    """
    match = re.findall(r"[a-zA-Z]\d\d", fname)
    match = [re.sub("0+(?!$)", "", loc) for loc in match]
    group_lst = ["T1"]
    split_fname = fname.split("_")
    group_quant = 0

    def screen_cond() -> tuple[list[str], int]:
        group_quant = 0

        if len(match) == 2:
            row = (match[0][0], match[1][0])
            col = (int(match[0][1:]), int(match[1][1:]))
            group_label = ""

            if row[0] == row[1]:
                group_quant = col[1] - col[0] + 2
                for i in range(group_quant - 1):
                    group_label = row[0] + str(col[0] + i)
                    group_lst.append(group_label)

            elif row[0] != row[1]:
                group_quant = 12 - col[0] + col[1] + 2
                for i in range(group_quant - 1):
                    if col[0] + i > 12:
                        group_label = row[1] + str(col[0] + i - 12)
                    elif col[0] + i <= 12:
                        group_label = row[0] + str(col[0] + i)
                    group_lst.append(group_label)

            else:
                raise ValueError(f"Error processing file: {fname}")

        # If the samples are not continues in the file, i.e., A01, A04, A09, B03
        elif len(match) > 2:
            group_quant = len(match) + 1
            group_lst.extend(match)

        elif len(match) < 2:
            raise ValueError(f"Check test and file name - {fname}")

        return group_lst, group_quant

    # Generate the groups for screening test type
    if split_fname[2] == "Screen":
        group_lst, group_quant = screen_cond()

    elif split_fname[2] == "Conc":
        if option == "under":
            conc_lst = [100, 75, 50, 25]
        elif option == "over":
            conc_lst = [100, 125, 150, 175, 200]
        else:
            raise ValueError(f"Unknown option - {option}")

        if match[0] == "A0":
            group_lst = ["T1_" + str(conc) for conc in conc_lst]
        else:
            name_lst = [
                f"{split_fname[2]}_{match[0]}_" + str(conc) for conc in conc_lst
            ]
            group_lst.extend(name_lst)
        group_quant = len(group_lst)

    elif split_fname[2] == "Ratio":
        trig = f"{split_fname[3]}_{match[0]}_{split_fname[5]}_{match[1]}"
        single_conc_lst = [100, 75, 50, 25, 0]
        ratio_tuple_lst = list(zip(single_conc_lst, single_conc_lst[::-1]))

        def combine(tuple_pair: tuple):
            return trig + "_" + str(tuple_pair[0]) + "_" + str(tuple_pair[-1])

        group_lst = list(map(combine, ratio_tuple_lst))
        group_quant = len(group_lst)

    return group_quant, group_lst


def group_generation(
    df: pd.DataFrame, con_tube_number: int
) -> tuple[list[list[str]], list[list]]:
    """
    Use to generate the groups based on the location of the control tube. Generate
    both general groups and presented groups. General groups are the groups that should be
    in the test file. The presented groups are the groups that are actually in the test
    file. This occurs when the data processing procedure observes negative or stagnet
    data points.

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe from the test file after the data processing.

    con_tube_number: int
        The location of the control tube in the test.

    Returns
    -------
    sample_group : list
        The list of test groups that should be in the test file.

    presented_group : list
        The list of test groups that are presented in the test file.

    """
    sorted_non_duplicate_tube_number = sorted(
        [int(i) for i in list(set(list(df.columns)))]
    )
    presented_tube_number = [str(i) for i in sorted_non_duplicate_tube_number]
    # Create subgroups of individual conditions (4 tubes per condition continuously)
    counter1 = 1
    counter2 = 1
    counter3 = 0
    number_of_samples = list(range(1, con_tube_number))
    sample_group = []
    while counter1 <= (con_tube_number - 1) / 4:
        subgroup = list(str())
        while counter3 != 4:
            subgroup.append(str(number_of_samples[counter2 - counter1]))
            counter2 += 1
            counter3 += 1
        sample_group.append(subgroup)
        counter1 += 1
        counter2 += 1
        counter3 = 0

    presented_groups = []
    for i in sample_group:
        presented_subgroups = []
        for j in i:
            if j in presented_tube_number:
                presented_subgroups.append(j)
        presented_groups.append(presented_subgroups)

    return sample_group, presented_groups


# TODO: Simplify this function
def data_combination(
    fn: str,
    ext: str,
    groups: int,
    cdata_path: str,
    opt: str | None = None,
) -> pd.DataFrame:
    """
    Combine the data across differnt trials with same experiemtnal
    condition into one dataframe

    Parameters
    ----------
    fn : str
        The based file name used for searching with glob.glob

    ext : str
        The file extension used for searching with glob.glob

    groups: intzr
        The number of groups in the trial

    cdata_path : str
        The save path for the combined dataframe in pickle format

    opt : str | None
        Optional augument with default None. Used in ratio study to indicate the
        direction of the ratio.

    Returns
    -------
    combined_df : pd.DataFrame
        The combined dataframe for the trials with the same testing conditions.

    """
    # initial data after subtraction of baseline and data cleaning
    count = 0
    print(fn)
    test_type = fn.split("_")[2]
    query_str = fn + f"*.{ext}"
    if test_type == "Conc":
        query_str = f"{fn}_{opt}_*.{ext}"

    combined_df = pd.DataFrame()
    time_lst_min = time_list_generation(60)
    for file in glob.glob(query_str):
        if test_type == "Screen":
            # Check if the file has the positive and negative control tube location reversed
            if file.split(".")[0].split("_")[-1] == "RC":
                pos_con_tube_number = (groups * 4) + 1
                neg_con_tube_number = pos_con_tube_number + 4
                col_read = [str(i) for i in range(1, neg_con_tube_number + 4)]
                reverse_control = True
                print(
                    f"This file contained the study with the reversed control order: {file}"
                )

            else:
                neg_con_tube_number = (groups * 4) + 1
                pos_con_tube_number = neg_con_tube_number + 4
                col_read = [str(i) for i in range(1, pos_con_tube_number + 4)]
                reverse_control = False

        elif test_type in {"Conc", "Ratio"}:
            neg_con_tube_number = (groups * 4) + 1
            pos_con_tube_number = neg_con_tube_number + 4
            col_read = [str(i) for i in range(1, pos_con_tube_number + 4)]
            reverse_control = False

        else:
            print(f"The current file name is {fn}")
            raise ValueError(f"The test type is unknown - {test_type}")

        col_read.insert(0, "Cycle")

        # Read the current csv file with the predefined column names
        print(f"currently reading file: {file}")
        df = pd.read_csv(file, usecols=col_read)
        df.insert(1, "time (min)", time_lst_min, True)
        df.set_index("Cycle", inplace=True)
        avg_neg_control = df.iloc[
            :, list(range(neg_con_tube_number, neg_con_tube_number + 4))
        ].mean(axis=1)
        pos_control = df.iloc[
            :, list(range(pos_con_tube_number, pos_con_tube_number + 4))
        ]
        pos_control_minus_baseline = pos_control.subtract(avg_neg_control, axis=0)
        avg_pos_control_minus_baseline = pos_control_minus_baseline.mean(axis=1)

        if reverse_control:
            data = df.iloc[:, list(range(1, pos_con_tube_number))]
        else:
            data = df.iloc[:, list(range(1, neg_con_tube_number))]

        data_minus_baseline = data.subtract(avg_neg_control, axis=0)

        # Finding the files that matches the condition,
        # combine them together and divide by the avg positive control

        norm_data = data_minus_baseline.div(avg_pos_control_minus_baseline, axis=0)

        final_value = norm_data.iloc[-1, :]
        ind_norm_data = norm_data.div(final_value, axis=1)

        # Drop the samples that did not started by cycle 3 and change in
        # values between initial and cycle 4 is less than 0.05

        drop_col = []
        for col in range(len(ind_norm_data.columns)):
            loc_0 = ind_norm_data.iloc[0, col]
            loc_3 = ind_norm_data.iloc[3, col]
            loc_4 = ind_norm_data.iloc[4, col]

            if (
                isinstance(loc_0, float)
                and isinstance(loc_3, float)
                and isinstance(loc_4, float)
            ):
                if loc_4 - loc_0 < 0.05 and loc_3 < 0:
                    drop_col.append(str(col + 1))
            else:
                raise ValueError()

        print(drop_col)
        norm_data.drop(columns=drop_col, inplace=True)

        # check if there are existing file for the condition.
        # If there is, append the new data to it, else create a new dataframe and export afterward.
        if count == 0:
            if not os.path.exists(cdata_path):
                combined_df = norm_data
            else:
                combined_df = pd.read_pickle(cdata_path)
                if not combined_df.equal(norm_data):
                    combined_df = pd.concat([combined_df, norm_data], axis=1)
        else:
            combined_df = pd.concat([combined_df, norm_data], axis=1)
        count += 1

        combined_df.to_pickle(cdata_path)

    return combined_df


def data_normalization(
    df: pd.DataFrame, groups: list[list], ndata_path: str
) -> pd.DataFrame:
    """
    Normalized data for each test file.

    Normalized the values from combined data to the standard release
    values and export it into csv and/or pickle file for storage and
    quicker access in Python.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe that combined all the same testing condition.

    groups: list
        The list of test groups that should be in the test file.

    ndata_path : str
        The save path for the dataframe in pickle format

    Returns
    -------
    standard_norm_df : pd.DataFrame
        The normalized values for each individual run in the test.

    """

    standard_value = df.loc[:, groups[0]]
    standard_value_average_max: float = standard_value.mean(axis=1).max()
    standard_norm_df = df / standard_value_average_max
    standard_norm_df.to_pickle(ndata_path)

    return standard_norm_df


def data_average(
    df: pd.DataFrame,
    time_list: list[float],
    group_dict: dict,
    presented_groups: list,
    sdata_path: str,
) -> pd.DataFrame:
    """
    Calculate basic stats for the normalized data and combined them into one dataframe

    Parameters
    ----------
    df : pd.DataFrame
        The normalized dataframe

    time_list : list[float]
        The list of time points correspondes to the individual measurment in minutes

    group_dict : dict
        The

    presented_groups : list
        The list of groups that are presented in the trials

    sdata_path : str
        The path for saving the combined and normalized dataframe in pickle format

    Returns
    -------
    norm_combined_data_df : pd.DataFrame
        The combined and normalized dataframe

    """

    norm_combined_data_df = pd.DataFrame()
    counter = 1

    for group in presented_groups:
        temp_col_name = f"{group_dict[counter]}"
        print(len(df.loc[:, group].columns))
        mean_col_name = temp_col_name + "_mean"
        std_col_name = temp_col_name + "_std"
        mean_col = df.loc[:, group].mean(axis=1)
        std_col = df.loc[:, group].std(axis=1)
        mean_col.rename(mean_col_name, inplace=True)
        std_col.rename(std_col_name, inplace=True)
        group_data = pd.concat([mean_col, std_col], axis=1)

        if counter == 1:
            norm_combined_data_df = group_data
        else:
            norm_combined_data_df = pd.concat(
                [norm_combined_data_df, group_data], axis=1
            )
        counter += 1

    norm_combined_data_df.insert(0, "time (min)", time_list)
    # export the combined file for storage + quick access
    norm_combined_data_df.to_pickle(sdata_path)

    return norm_combined_data_df


def data_summerization(path: str) -> None:
    """
    Combine the same test condition from screening study into one file.

    Parameters
    ----------
    path : str
        The location where the normalized pickle files are stored

    """
    os.chdir(path)
    t1_lst = []
    for fname in glob.glob("*" + "Screen" + "*_normalized.pkl"):
        print(fname)
        plate_num = fname.split("_")[3]
        groups, group_lst = group_query(fname)
        non_t1_cond = [f"{plate_num}_{i}" for i in group_lst if i != "T1"]

        df = pd.read_pickle(fname)
        _, presented = group_generation(df, (groups * 4 + 1))
        presented_dict = dict(zip(non_t1_cond, presented[1:]))

        for k, v in presented_dict.items():
            filtered_df = df.loc[:, v]
            filtered_df.to_pickle(f"./Individual/4WJ_HEX_Screen_{k}.pkl")

        t1_lst.append(df.loc[:, ["1", "2", "3", "4"]])

    t1_df = pd.concat(t1_lst, axis=1, ignore_index=False)
    t1_df.to_pickle("./Individual/4WJ_HEX_Screen_P0_T1.pkl")
    sum_t1_df = pd.DataFrame(
        dict(zip(["mean", "std"], [t1_df.mean(axis=1), t1_df.std(axis=1)]))
    )
    sum_t1_df.to_pickle("./Individual/4WJ_HEX_Screen_P0_T1_summarized.pkl")


def parameter(path: str) -> None:
    """
    Combine the fitted results into one DataFrame.

    Parameters
    ----------
    path : str
        The location where the fitted results are stored.

    """
    os.chdir(path)
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
