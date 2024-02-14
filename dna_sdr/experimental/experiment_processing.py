"""
This module is used to process the data from different test type

This module utilized the data_processing module to process the data collected for the different  
test type.

"""
import glob
import os
import shutil

import dna_sdr.experimental.data_processing as dp


# * Potentially add the option arg here for conc. study
def file_search(
    test_type: str, fluorophore: str, ext: str = "csv"
) -> tuple[dict[int, str], list[str]]:
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
    fname_search = f"4WJ_{fluorophore}_{test_type}_*.{ext}"
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
    test_type: str, trigger_selected: str, ext: str, fluorophore: str
) -> None:
    """
    xxx

    Parameters
    ----------
    test_type : str

    trigger_selected : str

    """

    # Save file path generation
    # paths[0] = fname
    # paths[1] = combined_data_path
    # paths[2] = norm_cdata_path
    # paths[3] = sum_data_path
    # paths[4] = processed_path
    paths = dp.file_path_generation(test_type, trigger_selected, fluorophore)
    time_list_mins = dp.time_list_generation(60)

    if test_type not in {"Conc", "Screen", "Ratio"}:
        raise ValueError()

    if test_type == "Conc":
        option = trigger_selected.split("_")[-1]
        groups, group_list = dp.group_query(paths[0], option)
    else:
        groups, group_list = dp.group_query(paths[0])

    # Combined the data together and check for conditions that are not accurate
    # after subtracting the negative control and generate a combined data dataframe for
    # further analysis.
    combined_data_df = dp.data_combination(paths[0], ext, groups, paths[1])

    # Generate the groups presented in the combined data dataframe
    group_dict = dict(zip(list(range(1, len(group_list) + 1)), group_list))
    _, presented_groups = dp.group_generation(combined_data_df, (groups * 4 + 1))

    # Normalized the dataframe fo combined data with standard release value (T1) and export
    # into a pickle file for storage.
    norm_data_df = dp.data_normalization(combined_data_df, presented_groups, paths[2])

    # Combine the values from norm_data_df to calculate average and standard deviation of
    # each condition and combined it all into one dataframe.
    dp.data_average(
        norm_data_df, time_list_mins, group_dict, presented_groups, paths[3]
    )

    # Move the processed data to the processed folder in specific folder
    # corresponding to the test + conditions if the folder does not exist,
    # create the folder and move the file there, otherwise, move the file to the
    # corresponding folder
    if not os.path.exists(paths[4]):
        os.mkdir(paths[4])
    else:
        print(f"Folder name: {paths[4]} already exists.")

    for f in glob.glob(paths[0] + f"*.{ext}"):
        shutil.move(os.path.join(os.getcwd(), f), os.path.join(paths[4], f))


def main(test: str, fluorophore: str):
    """
    xxx

    Parameters
    ----------
    test : str


    """
    t_dict, _ = file_search(test, fluorophore)
    for trig in t_dict.values():
        print(trig)
        exp_data_processing(test, trig, "csv", fluorophore)


if __name__ == "__main__":
    TEST = "Ratio"
    F = "HEX"
    main(TEST, F)
