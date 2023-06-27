import glob
import os
import pandas as pd


# TODO: Update the data combination function for concentration screening
def data_combination(fn, ext, groups, t_list_mins, cdata_path, opt=None):
    # initial data after subtraction of baseline and data cleaning
    count = 0
    combined_df = pd.DataFrame
    print(fn)
    test_type = fn.split("_")[2]
    query_str = fn + "*.{}".format(ext)
    if test_type == "Conc":
        query_str = "{}_{}_*.{}".format(fn, opt, ext)

    for file in glob.glob(query_str):
        if test_type == "Screen":
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

        elif test_type == "Conc":
            neg_con_tube_number = (groups * 4) + 1
            pos_con_tube_number = neg_con_tube_number + 4
            col_read = [str(i) for i in range(1, pos_con_tube_number + 4)]
            reverse_control = False

        else:
            print("The current file name is {}".format(fn))
            raise ValueError("The test type is unknown - {}".format(test_type))

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
                combined_df = norm_data
            else:
                combined_df = pd.read_pickle(cdata_path)
                if not combined_df.equal(norm_data):
                    combined_df = pd.concat([combined_df, norm_data], axis=1)
        else:
            combined_df = pd.concat([combined_df, norm_data], axis=1)
        count += 1

    return combined_df


# TODO: Check if it works for the concentration study prior to uncomment pickle export
def data_normalization(df, groups, ndata_path):
    """
    Normalized the values from combined data to the standard release
    values and export it into csv and/or pickle file for storage and
    quicker access in Python
    """
    standard_value = df.loc[:, groups[0]]
    standard_value_average_max: float = standard_value.mean(axis=1).max()
    standard_norm_df = df / standard_value_average_max
    standard_norm_df.to_pickle(ndata_path)

    return standard_norm_df


# TODO: Check if it works for the concentration study prior to uncomment pickle export
def data_average(df, time_list, group_dict, presented_groups, sdata_path):
    # Calculate the average and stdev of the normalized combined data
    # and combined them into one dataframe
    norm_combined_data_df = pd.DataFrame
    counter = 1

    for group in presented_groups:
        temp_col_name = "{}".format(group_dict[counter])
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

    norm_combined_data_df.insert(0, "time (min)", time_list, True)
    print(norm_combined_data_df)

    # export the combined file for storage + quick access
    norm_combined_data_df.to_pickle(sdata_path)

    return
