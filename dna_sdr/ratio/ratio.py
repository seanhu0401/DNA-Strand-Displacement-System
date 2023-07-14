import glob
import os
import shutil
import regex as re
import dna_sdr.data_process.list_generation as lst_gen
import dna_sdr.data_process.group as grouping
import dna_sdr.data_process.data_processing as dp

if __name__ == "__main__":
    cwd = os.getcwd()
    path = cwd + "/dna_sdr/IO/Data_Input/"
    extension = "csv"
    os.chdir(path)

    test_type = "Ratio"
    file_list = []

    # Loop through the files from the input folder and find the matching files for processing
    # Add the trigger type into a list and convert to a dict for selection
    fname_search = "4WJ_HEX_" + "{}_*.{}".format(test_type, extension)
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
        path[0] = partial fname
        path[1] = combined_data_path
        path[2] = norm_cdata_path
        path[3] = sum_data_path
        path[4] = processed_path
        """
        paths = lst_gen.file_path_generation(test_type, trig_type_selected)
        fn = [f for f in file_list if re.match(paths[0] + "_*", f)]

        groups, group_list = grouping.group_query(paths[0])
        group_dict = dict(zip(list(range(1, len(group_list) + 1)), group_list))
        time_list_mins = lst_gen.time_list_generation(60)
