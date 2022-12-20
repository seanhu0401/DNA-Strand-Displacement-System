import os
import glob
import pandas as pd
import regex as re


def trig_list_gen(fname_list):
    trig_list = []
    re_query = "[a-zA-Z]\d\d"
    for fname in fname_list:
        split_name = fname.split(".")[0].split("_")
        plate_number = split_name[3]

        match = re.findall(re_query, fname)
        trig_type = plate_number + "_" + "_".join(match)

        if trig_type not in trig_list:
            trig_list.append(trig_type)

    return trig_list


def file_path_generation(t_type, trig_type):
    pickle_extension = "pkl"
    pickle_save_path = "../Output/Pickles"
    parent_dir = "../Processed/"

    if t_type == "Ratio" or t_type == "Screen" or t_type == "Conc":
        dir_name = "4WJ_HEX_{t_test}_{Trig}".format(t_test=t_type, Trig=trig_type)

    else:
        raise Exception("unknown test type")

    combined_data_after_baseline_subtraction_pickle = "{base}.{extension}".format(
        base=dir_name, extension=pickle_extension
    )
    norm_combined_data_pickle = "{base}_normalized.{extension}".format(
        base=dir_name, extension=pickle_extension
    )
    summerized_data_pickle = "{base}_summerized.{extension}".format(
        base=dir_name, extension=pickle_extension
    )
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


def group_query(fname):
    re_query = "[a-zA-Z]\d\d"
    match = re.findall(re_query, fname)
    match = [re.sub("0+(?!$)", "", loc) for loc in match]
    group_lst = ["T1"]

    # If the samples are contiunes in the file,i.e., A01-A03 or A09-B03
    if len(match) <= 2:
        row = (match[0][0], match[1][0])
        col = (int(match[0][1:]), int(match[1][1:]))

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
            raise Exception("Error processing file: {}".format(fname))

    # If the samples are not continues in the file, i.e., A01, A04, A09, B03
    elif len(match) > 2:
        group_quant = len(match) + 1
        group_lst.extend(match)

    return group_quant, group_lst


def time_list_generation(num_cycle, init_time_sec=231):
    cycle_list = list(range(1, num_cycle))
    t_list_seconds = [init_time_sec]

    for cycle in cycle_list:
        t = t_list_seconds[-1] + 20 + 5 * cycle
        t_list_seconds.append(t)

    t_list_mins = [round(t / 60, 2) for t in t_list_seconds]

    return t_list_mins


def data_combination(fn, ext, groups, t_list_mins, cdata_path):
    # initial data after subtraction of baseline and data cleaning
    count = 0
    cdata = pd.DataFrame
    print(fn)

    for file in glob.glob(fn + "*.{}".format(ext)):
        # Check if the file has the positive and negative control tube location reversed
        if file.split(".")[0].split("_")[-1] == "RC":
            pos_con_tube_number = (groups * 4) + 1
            neg_con_tube_number = pos_con_tube_number + 4
            col_read = [str(i) for i in range(1, neg_con_tube_number + 4)]
            reverse_control = True
            print(
                "This file contained the study with the reversed control order: {}".format(
                    file
                )
            )

        else:
            neg_con_tube_number = (groups * 4) + 1
            pos_con_tube_number = neg_con_tube_number + 4
            col_read = [str(i) for i in range(1, pos_con_tube_number + 4)]
            reverse_control = False

        col_read.insert(0, "Cycle")

        # Read the current csv file with the predefined column names
        print("currently reading file: {}".format(file))
        df = pd.read_csv(file, usecols=col_read)
        df.insert(1, "time (min)", t_list_mins, True)
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

        # Finding the files that matches the condition, combine them together and divide by the avg positive control
        norm_data = data_minus_baseline.div(avg_pos_control_minus_baseline, axis=0)

        final_value = norm_data.iloc[-1, :]
        ind_norm_data = norm_data.div(final_value, axis=1)

        # Drop the samples that did not started by cycle 3 and change in values between initial and cycle 4 is less than 0.05
        drop_col = list()
        for col in range(len(ind_norm_data.columns)):
            value_change = ind_norm_data.iloc[4, col] - ind_norm_data.iloc[0, col]
            if value_change < 0.05 and ind_norm_data.iloc[3, col] < 0:
                drop_col.append(str(col + 1))
        print(drop_col)
        norm_data.drop(columns=drop_col, inplace=True)

        # check if there are existing file for the condition. If there is, append the new data to it, else create a new dataframe and export afterward.
        if count == 0:
            if not os.path.exists(cdata_path):
                cdata = norm_data
            else:
                cdata = pd.read_pickle(cdata_path)
                cdata = pd.concat([cdata, norm_data], axis=1)
        else:
            cdata = pd.concat([cdata, norm_data], axis=1)
        count += 1

    return cdata


def group_generation(df, con_tube_number):
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
        subgroup = []
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
