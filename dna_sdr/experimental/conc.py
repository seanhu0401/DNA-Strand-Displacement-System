import glob
import os
import shutil
import regex as re
import dna_sdr.experimental.data_processing as dp

if __name__ == "__main__":
    cwd = os.getcwd()
    path = cwd + "/dna_sdr/IO/Data_Input/"
    extension = "csv"
    os.chdir(path)

    test_type = "Conc"
    file_list = []

    # Loop through the files from the input folder and find the matching files for processing
    # Add the trigger type into a list and convert to a dict for selection
    fname_search = "4WJ_HEX_" + f"{test_type}_*.{extension}"
    for f in glob.glob(fname_search):
        if f not in file_list:
            file_list.append(f)

    trig_list = dp.trig_list_gen(file_list)
    trig_type_dict = dict(zip(list(range(1, len(trig_list) + 1)), trig_list))

    if not trig_list:
        print("No file found for processing")

    else:
        print(trig_type_dict)
        trig_type_input = int(input("Type of trigger:") or "1")
        trig_type_selected = trig_type_dict.get(trig_type_input)

        """
        Save file path generation
        path[0] = partial fname
        path[1] = combined_data_path
        path[2] = norm_cdata_path
        path[3] = sum_data_path
        path[4] = processed_path
        """
        paths = dp.file_path_generation(test_type, trig_type_selected)
        fn = [f for f in file_list if re.match(paths[0] + "_*", f)]

        lower_len = len([f for f in fn if re.match(r"\w*_under", f)])
        upper_len = len([f for f in fn if re.match(r"\w*_over", f)])
        if len(fn) == lower_len:
            option = "under"
        elif len(fn) == upper_len:
            option = "over"
        groups, group_list = dp.group_query(paths[0], option)
        group_dict = dict(zip(list(range(1, len(group_list) + 1)), group_list))
        time_list_mins = dp.time_list_generation(60)

        """
        Combined the data together and check for conditions that are not accurate
        after normalization with the positive and negative control and generate a
        combined data dataframe for further analysis.
        """

        combined_data_df = dp.data_combination(
            paths[0], extension, groups, paths[1], opt=option
        )

        # export the combined file for storage + quick access
        combined_data_df.to_pickle(path[1])

        # Generate the groups presented in the combined data dataframe
        con_tube_number = groups * 4 + 1
        group_of_samples, presented_groups = dp.group_generation(
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
        if not os.path.exists(paths[4]):
            os.mkdir(paths[4])
        else:
            print("Folder name: {} already exists.".format(paths[4]))

        for f in glob.glob(paths[0] + "*.{}".format(extension)):
            shutil.move(path + f, paths[4] + "/" + f)
