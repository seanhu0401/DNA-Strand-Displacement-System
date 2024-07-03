"""
Module docstring
"""

from copy import deepcopy
import os

import pandas as pd
import regex as re
from more_itertools import consecutive_groups

from dna_sdr.sequence.dna_utilits import DNA, Trigger

comp_dict = {"t": 1, "g": 2, "c": 3, "a": 4}
base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}
# TODO: Update the mismatch conversion dict
mismatch_type_dict = {
    "a-t": 1,
    "a-g": 2,
    "a-c": 3,
    "t-a": 4,
    "t-c": 5,
    "t-g": 6,
    "c-a": 7,
    "c-t": 8,
    "c-g": 9,
    "g-t": 10,
    "g-c": 11,
    "g-a": 12,
}


def find_range(iterable: list[int]):
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    for group in consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0]
        else:
            yield group[0], group[-1]


# *Think of better name for the method
def sequence_comparison(seq_1: str, seq_2: str) -> dict[int, int]:
    """
    This function is currently used to compare to sequence, mainly the triggers, to identify the
    location(s) that has/have mismatch(es). They should be equal length in bases if there is no
    modification to the toehold + dangles - not applicable with those modification currently

    Parameters
    ----------

    Returns
    -------

    """
    comparison_lst: list[int] = []
    base_loc_lst: list[int] = []
    counter: int = 0

    try:
        strand_1 = DNA(seq_1)
        strand_2 = DNA(seq_2)
    except ValueError as e:
        raise TypeError("The input is not all valid nucleotide") from e

    if strand_1.base_count != strand_2.base_count:
        raise ValueError("The two sequences are not equal length")

    while counter < strand_1.base_count:
        comp_val: int = 0
        if strand_1.seq[counter] == strand_2.seq[counter]:
            comparison_lst.append(comp_val)
        elif strand_1.seq[counter] != strand_2.seq[counter]:
            comp_val = comp_dict[strand_2.seq[counter]]
            comparison_lst.append(comp_val)
        base_loc_lst.append(counter + 1)
        counter += 1
    seq_comp_dict = dict(zip(base_loc_lst, comparison_lst))

    return seq_comp_dict


def trig_alignment(
    trigger1: str, trigger2: str, toehold_b: int = 7, overhang_b: int = 5
) -> tuple[int, int, int, list[int], list[int]]:
    """
    xxx

    trigger1 is the complementory trigger for the system
    trigger2 is the modified trigger for the system
    input toehold and overhang bases if the inital trigger does not have
    7 toehold bases and 5 overhang bases.

    Parameters
    ----------


    Returns
    -------


    """

    def tb_ob_calculation(
        t1: DNA, t2: DNA, toehold_b: int, overhang_b: int
    ) -> tuple[int, int]:
        base_difference = abs(t1.base_count - t2.base_count)
        new_tb = toehold_b
        new_ob = overhang_b
        if toehold_b == 7 and overhang_b == 5:
            if base_difference < overhang_b and base_difference != 0:
                new_ob: int = overhang_b - base_difference
            elif (toehold_b + overhang_b) >= base_difference >= overhang_b:
                new_tb: int = toehold_b + overhang_b - base_difference
                new_ob: int = 0
            elif base_difference > (toehold_b + overhang_b):
                raise ValueError(
                    "Missing toehold region and overhang - please check input sequence"
                )
        return new_tb, new_ob

    try:
        trig1 = DNA(trigger1)
        trig2 = DNA(trigger2)
    except ValueError as e:
        raise TypeError("Not all sequence given is DNA sequences") from e

    counter = tb_mismatch = ob_mismatch = mismatch = 0
    mismatch_loc: list[int] = []
    mismatch_type: list[int] = []

    new_tb, new_ob = tb_ob_calculation(trig1, trig2, toehold_b, overhang_b)

    for base in trig2.reverse_seq():
        base_at_loc = trig1.reverse_seq()[counter]
        if base != base_at_loc and counter < (trig2.base_count - new_ob):
            location = trig1.base_count - counter
            mismatch_loc.append(location)
            mismatch_str = f"{base_at_loc}-{base}"
            mismatch_type.append(mismatch_type_dict[mismatch_str])
            mismatch += 1
            if location <= (new_tb + new_ob):
                if location > new_ob:
                    tb_mismatch += 1
                else:
                    ob_mismatch += 1
        counter += 1

    new_tb = new_tb - tb_mismatch
    new_ob = new_ob - ob_mismatch
    mismatch_loc.sort()
    mismatch_type.reverse()

    return new_tb, new_ob, mismatch, mismatch_loc, mismatch_type


def name_generation(
    toehold_b: int, overhang_b: int, mismatch: list[int], type_mismatch: list[int]
) -> str:
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    mismatch_groups = list(find_range(mismatch))

    base_name = f"{toehold_b}tb_5'{overhang_b}ob"

    full_mismatch_name = str()

    for loc, group in enumerate(mismatch_groups):
        if isinstance(group, int):
            base_number = 1
            mismatch_name = f"_{group}_{base_number}mb_{type_mismatch[loc]}"
        else:
            base_number = group[-1] - group[0] + 1
            type_group = tuple(type_mismatch[:base_number])
            mismatch_name = f"_{group[0]}_{base_number}mb_{type_group}"
            del type_mismatch[:base_number]

        full_mismatch_name = full_mismatch_name + mismatch_name

    return base_name + full_mismatch_name


def seq_info(fname: str) -> pd.DataFrame:
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    xls = pd.ExcelFile(fname)
    master_list = []

    for sname in xls.sheet_names:
        if sname == "T3":
            trig = "CAC CAC CCA TCT CAA AAC TCA TCT CAT CCA ACA"
        else:
            trig = "CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"
        dna_seq_df = pd.read_excel(fname, sheet_name=sname)
        loc_df = dna_seq_df.filter(regex=r"P\d Location", axis=1)
        dna_seq = dna_seq_df["Seqences (5' to 3')"]
        for count in range(len(dna_seq)):
            location = loc_df.iloc[count].dropna()
            for well in location:
                plate = re.findall(r"P\d", location.index[0])[0]
                plate_loc = "_".join([plate, well])
            seq = dna_seq.iloc[count]
            info = trig_alignment(trig, seq)
            info_copy = deepcopy(info)
            name = name_generation(*info_copy[:2], *info_copy[3:])
            if sname == "T3":
                name = "_".join(["T3", name])
            try:
                dna_trig = Trigger(seq, name, *info, plate_loc)
            except UnboundLocalError:
                dna_trig = Trigger(seq, name, *info)

            master_list.append(dna_trig)

    trig_df = pd.DataFrame([vars(trig) for trig in master_list])

    return trig_df


if __name__ == "__main__":
    FNAME = "./dna_sdr/DNA_Strands.xlsx"
    PICKLE_NAME = "./dna_sdr/pickles/trig_info.pkl"
    info_df = seq_info(FNAME)

    max_mismatch: int = max(info_df["mismatch"])
    location_str_lst = []
    type_str_lst = []
    for count in range(max_mismatch):
        location_str_lst.append(f"mismatch_location_{count+1}")
        type_str_lst.append(f"mismatch_type_{count+1}")
    info_df[location_str_lst] = pd.DataFrame(
        info_df.mismatch_loc.to_list(), index=info_df.index
    )
    info_df[type_str_lst] = pd.DataFrame(
        info_df.mismatch_type.to_list(), index=info_df.index
    )
    print(info_df[["name", "plate_loc"]])

    if os.path.exists(PICKLE_NAME):
        trigger_df = pd.read_pickle(PICKLE_NAME)
        if info_df.equals(trigger_df):
            print("Same dataframe - no new file")
        else:
            info_df.to_pickle(PICKLE_NAME)
    else:
        info_df.to_pickle(PICKLE_NAME)
