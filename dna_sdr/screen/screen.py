import glob
import os
import dna_sdr.data_process.list_generation as lst_gen
import dna_sdr.data_process.data_processing as dp
import dna_sdr.data_process.group as grouping
import shutil

if __name__ == "__main__":
    cwd = os.getcwd()
    path = cwd + "/dna_sdr/IO/Data_Input/"
    extension = "csv"
    os.chdir(path)

    test_type = "Screen"
    file_list = []

    # Loop through the files from the input folder and find the matching files for processing
    # Add the trigger type into a list and convert to a dict for selection
    fname_search = "4WJ_HEX_Screen" + "*.{}".format(extension)
    for f in glob.glob(fname_search):
        if f not in file_list:
            file_list.append(f)

    trig_list = lst_gen.trig_list_gen(file_list)
    trig_type_dict = dict(zip(list(range(1, len(trig_list) + 1)), trig_list))

    if not trig_list:
        print("No file found for processing")

    else:
        print(trig_type_dict)
        trig_type_input = int(input("Type of trigger:") or "1")
        trig_type_selected = trig_type_dict.get(trig_type_input)

        """
        Save file path generation
        path[0] = fname 
        path[1] = combined_data_path
        path[2] = norm_cdata_path
        path[3] = sum_data_path
        path[4] = processed_path
        """
        paths = lst_gen.file_path_generation(test_type, trig_type_selected)

        groups, group_list = grouping.group_query(paths[0])
        group_dict = dict(zip(list(range(1, len(group_list) + 1)), group_list))
        time_list_mins = lst_gen.time_list_generation(60)

        """
        Combined the data together and check for conditions that are not accurate
        after normalization with the positive and negative control and generate a
        combined data dataframe for further analysis.
        """
        combined_data_df = dp.data_combination(
            path[0], extension, groups, time_list_mins, path[1]
        )

        # export the combined file for storage + quick access
        combined_data_df.to_pickle(path[1])

        # Generate the groups presented in the combined data dataframe
        con_tube_number = groups * 4 + 1
        group_of_samples, presented_groups = grouping.group_generation(
            combined_data_df, con_tube_number
        )

        """
        Normalized the dataframe fo combined data with standard release value (T1) and 
        export the dataframe into csv and/or pickle files for storage and quicker access
        in Python
        """
        norm_data_df = dp.data_normalization(
            combined_data_df, group_of_samples, paths[2]
        )

        """
        Combine the values from norm_data_df to calculate average and standard devaition of
        each condition and combined it all into one dataframe. Export the dataframe into 
        csv and/or pickle files for storage and quicker access in Python
        """
        dp.data_average(
            norm_data_df, time_list_mins, group_dict, presented_groups, paths[3]
        )

        """
        Move the processed data to the processed folder in specific folder 
        corresponding to the test + conditions if the folder does not exist, 
        create the folder and move the file there, otherwise, move the file to the 
        corresponsing folder
        """
        if not os.path.exists(path[4]):
            os.mkdir(path[4])
        else:
            print("Folder name: {} already exists.".format(path[4]))

        for f in glob.glob(path[0] + "*.{}".format(extension)):
            shutil.move(path + f, path[4] + "/" + f)


else:
    print()
