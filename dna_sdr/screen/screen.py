import glob
import os
import dna_sdr.data_process.data_organization as do
import pandas as pd
import shutil

if __name__ == "__main__":
    cwd = os.getcwd()
    path = cwd + "/dna_sdr/IO/Data_Input/"
    extension = "csv"
    os.chdir(path)

    test_type = "Screen"
    file_list = []

    # Loop through the files from the input folder and find the matching files for processing
    # Add the trigger type into a list and convert to a dict for manual selection
    fname_search = "4WJ_HEX_Screen" + "*.{}".format(extension)
    for f in glob.glob(fname_search):
        if f not in file_list:
            file_list.append(f)

    trig_list = do.trig_list_gen(file_list)
    trig_type_dict = dict(zip(list(range(1, len(trig_list) + 1)), trig_list))

    if not trig_list:
        print("No file found for processing")

    else:
        print(trig_type_dict)
        trig_type_input = int(input("Type of trigger:") or "1")
        trig_type_selected = trig_type_dict.get(trig_type_input)

        # Save file path generation
        (
            fname,
            combined_data_path,
            norm_cdata_path,
            sum_data_path,
            processed_path,
        ) = do.file_path_generation(test_type, trig_type_selected)

        groups, group_list = do.group_query(fname)
        group_dict = dict(zip(list(range(1, len(group_list) + 1)), group_list))
        time_list_mins = do.time_list_generation(60)

        # Combined the data together and check for conditions that are not accurate after normalization with the positive and negative control - generate a combined data dataframe for further analysis.
        combined_data = do.data_combination(
            fname,
            extension,
            groups,
            time_list_mins,
            combined_data_path,
        )

        # export the combined file for storage + quick access
        combined_data.to_pickle(combined_data_path)
        # Generate the groups presented in the combined data
        con_tube_number = groups * 4 + 1
        group_of_samples, presented_groups = do.group_generation(
            combined_data, con_tube_number
        )

        # Normalized the values from combined data to the T1 values and export it into csv and pickle files.
        T1_values = combined_data.loc[:, group_of_samples[0]]
        T1_value_average_max: float = T1_values.mean(axis=1).max()
        T1_norm_combined_data = combined_data / T1_value_average_max
        T1_norm_combined_data.to_pickle(norm_cdata_path)

        # Calculate the average and stdev of the normalized combined data and combined them into one dataframe
        norm_combined_data = pd.DataFrame
        counter = 1

        for group in presented_groups:
            temp_col_name = "{}".format(group_dict[counter])
            print(len(T1_norm_combined_data.loc[:, group].columns))
            mean_col_name = temp_col_name + "_mean"
            std_col_name = temp_col_name + "_std"
            mean_col = T1_norm_combined_data.loc[:, group].mean(axis=1)
            std_col = T1_norm_combined_data.loc[:, group].std(axis=1)
            mean_col.rename(mean_col_name, inplace=True)
            std_col.rename(std_col_name, inplace=True)
            group_data = pd.concat([mean_col, std_col], axis=1)

            if counter == 1:
                norm_combined_data = group_data
            else:
                norm_combined_data = pd.concat([norm_combined_data, group_data], axis=1)
            counter += 1

        norm_combined_data.insert(0, "time (min)", time_list_mins, True)

        # export the combined file for storage + quick access
        norm_combined_data.to_pickle(sum_data_path)

        # Move the processed data to the processed folder in specific folder corresponding to the test + conditions - if the folder does not exist, create the folder and move the file there, otherwise, move the file to the corresponsing folder
        if not os.path.exists(processed_path):
            os.mkdir(processed_path)
        else:
            print("Folder name: {} already exists.".format(processed_path))

        for f in glob.glob(fname + "*.{}".format(extension)):
            shutil.move(path + f, processed_path + "/" + f)


else:
    print()
