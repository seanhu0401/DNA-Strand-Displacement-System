import glob
import os
import pandas as pd


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
