"""
This module is used to process the data from differnet test type

This module utilized the data_processing module to process the data collected for the different  
test type.

"""
import glob
import os
import shutil
import regex as re
import pandas as pd
import dna_sdr.experimental.data_processing as dp


def file_search(ext: str, test_type: str) -> tuple[dict[int, str], list[str]]:
    """
    xxx

    Parameters
    ----------
    ext : str

    test_type : str

    Returns
    -------
    trig_type_dict : dict[int, str]

    trig_list : list[str]

    """
    cwd = os.getcwd()
    path = cwd + "/dna_sdr/IO/Data_Input/"
    os.chdir(path)

    file_list = []

    # Loop through the files from the input folder and find the matching files for processing
    # Add the trigger type into a list and convert to a dict for selection
    fname_search = "4WJ_HEX_" + f"{test_type}_*.{ext}"
    for f in glob.glob(fname_search):
        if f not in file_list:
            file_list.append(f)

    trig_list = dp.trig_list_gen(file_list)
    if not trig_list:
        raise ValueError("No file found for processing")

    trig_type_dict = dict(zip(list(range(1, len(trig_list) + 1)), trig_list))

    return trig_type_dict, trig_list


def trig_selection(trig_type_dict: dict[int, str]) -> str:
    """
    xxx

    Parameters
    ----------
    trig_type_dict : dict[int, str]

    Returns
    -------
    selected_trig : str

    """
    print(trig_type_dict)
    trig_type_input = int(input("Type of trigger:") or "1")
    selected_trig = trig_type_dict.get(trig_type_input)

    if not isinstance(selected_trig, str):
        raise ValueError

    return selected_trig


def exp_data_processing(
    test_type: str, trigger_selected: str, file_list: list[str], ext: str
) -> None:
    """
    xxx

    Parameters
    ----------
    test_type : str

    trigger_selected : str

    Returns
    -------
    selected_trig : str

    """

    def _concentration_study() -> tuple[int, list[str], pd.DataFrame]:
        fn_lst = [f for f in file_list if re.match(paths[0] + "_*", f)]

        lower_len = len([f for f in fn_lst if re.match(r"\w*_under", f)])
        upper_len = len([f for f in fn_lst if re.match(r"\w*_over", f)])
        if len(fn_lst) == lower_len:
            option = "under"
        elif len(fn_lst) == upper_len:
            option = "over"
        else:
            raise ValueError()
        groups, group_list = dp.group_query(paths[0], option)
        combined_data_df = dp.data_combination(
            paths[0], ext, groups, paths[1], opt=option
        )
        return groups, group_list, combined_data_df

    # Save file path generation
    # path[0] = fname
    # path[1] = combined_data_path
    # path[2] = norm_cdata_path
    # path[3] = sum_data_pathghfj
    # path[4] = processed_path
    paths = dp.file_path_generation(test_type, trigger_selected)
    time_list_mins = dp.time_list_generation(60)

    if not test_type in {"Conc", "Screen", "Ratio"}:
        raise ValueError()

    if test_type == "Conc":
        groups, group_list, combined_data_df = _concentration_study()
    else:
        groups, group_list = dp.group_query(paths[0])
        # Combined the data together and check for conditions that are not accurate
        # after normalization with the positive and negative control and generate a
        # combined data dataframe for further analysis.
        combined_data_df = dp.data_combination(paths[0], ext, groups, paths[1])

    group_dict = dict(zip(list(range(1, len(group_list) + 1)), group_list))

    # Generate the groups presented in the combined data dataframe
    group_of_samples, presented_groups = dp.group_generation(
        combined_data_df, (groups * 4 + 1)
    )

    # Normalized the dataframe fo combined data with standard release value (T1) and
    # export the dataframe into csv and/or pickle files for storage and quicker access
    # in Python
    norm_data_df = dp.data_normalization(combined_data_df, group_of_samples, paths[2])

    # Combine the values from norm_data_df to calculate average and standard devaition of
    # each condition and combined it all into one dataframe. Export the dataframe into
    # csv and/or pickle files for storage and quicker access in Python
    dp.data_average(
        norm_data_df, time_list_mins, group_dict, presented_groups, paths[3]
    )

    # Move the processed data to the processed folder in specific folder
    # corresponding to the test + conditions if the folder does not exist,
    # create the folder and move the file there, otherwise, move the file to the
    # corresponsing folder
    if not os.path.exists(paths[4]):
        os.mkdir(paths[4])
    else:
        print(f"Folder name: {paths[4]} already exists.")

    for f in glob.glob(paths[0] + f"*.{ext}"):
        shutil.move(os.getcwd() + f, paths[4] + "/" + f)
