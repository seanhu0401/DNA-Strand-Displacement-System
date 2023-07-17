import os
import regex as re


def trig_list_gen(fname_list):
    trig_list = []
    re_query = "[a-zA-Z]\d\d"
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


def file_path_generation(test_type, trig_type):
    pickle_extension = "pkl"
    pickle_save_path = "../Output/Pickles"
    parent_dir = "../Processed/"

    if test_type == "Ratio" or test_type == "Screen" or test_type == "Conc":
        dir_name = "4WJ_HEX_{t_test}_{Trig}".format(t_test=test_type, Trig=trig_type)

    else:
        raise Exception("unknown test type - {}".format(test_type))

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


def time_list_generation(num_cycle, init_time_sec=231):
    cycle_list = list(range(1, num_cycle))
    t_list_seconds = [init_time_sec]

    for cycle in cycle_list:
        t = t_list_seconds[-1] + 20 + 5 * cycle
        t_list_seconds.append(t)

    t_list_mins = [round(t / 60, 2) for t in t_list_seconds]

    return t_list_mins
