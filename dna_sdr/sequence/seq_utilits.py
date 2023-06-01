from dna_sdr.sequence.dna_utilits import DNA, trigger
import pandas as pd
import regex as re
import os
from more_itertools import consecutive_groups

comp_dict = {"t": 1, "g": 2, "c": 3, "a": 4}
base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}


def find_range(iterable):
    for group in consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0]
        else:
            yield group[0], group[-1]


# *Think of better name for the method
# This function is currently used to compare to sequence, mainly the triggers, to identify the location(s) that has/have mismatch(es)
# They should be equal length in bases if there is no modification to the toehold + dangles - not applicable with those modification currently
def sequence_comparison(seq_1: str, seq_2: str):
    comparison_lst = []
    base_loc_lst = []
    counter = 0

    try:
        strand_1 = DNA(seq_1)
        strand_2 = DNA(seq_2)
    except ValueError:
        raise TypeError("The input is not a DNA type")

    if strand_1.base_count != strand_2.base_count:
        raise Exception("The two sequences are not equal length")

    else:
        while counter < strand_1.base_count:
            comp_val = 0
            if strand_1.seq[counter] == strand_2.seq[counter]:
                comparison_lst.append(comp_val)
            elif strand_1.seq[counter] != strand_2.seq[counter]:
                comp_val = comp_dict[strand_2.seq[counter]]
                comparison_lst.append(comp_val)
            base_loc_lst.append(counter + 1)
            counter += 1
        seq_comp_dict = dict(zip(base_loc_lst, comparison_lst))

    return seq_comp_dict


def trig_aligment(trigger1: str, trigger2: str, tb=7, ob=5):
    """
    trigger1 is the complementory trigger for the system
    trigger2 is the modified trigger for the system
    input toehold and overhang bases if the inital trigger does not have
    7 toehold bases and 5 overhang bases.
    """
    try:
        trig1 = DNA(trigger1)
        trig2 = DNA(trigger2)
    except ValueError:
        raise TypeError("Not all sequence given is DNA sequences")

    counter = tb_mismatch = ob_mismatch = mismatch = 0
    mismatch_loc = list()
    mismatch_type = list()

    base_difference = abs(trig1.base_count - trig2.base_count)
    if tb == 7 and ob == 5:
        if base_difference < ob and base_difference != 0:
            new_ob = ob - base_difference
        elif (base_difference >= ob) and (base_difference <= (tb + ob)):
            new_tb = tb + ob - base_difference
            new_ob = 0
        elif base_difference == 0:
            new_tb = tb
            new_ob = ob
        elif base_difference > (tb + ob):
            raise ValueError(
                "Missing toehold region and overhang - please check input sequence"
            )
    else:
        new_tb = tb
        new_ob = ob

    for base in trig2.reverse_seq():
        base_at_loc = trig1.reverse_seq()[counter]

        if base == base_at_loc and counter < (trig2.base_count - new_ob):
            pass
        elif base != base_at_loc and counter < (trig2.base_count - new_ob):
            location = trig1.base_count - counter
            mismatch_loc.append(location)
            mismatch_type.append("{}-{}".format(base_at_loc, base))
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


def name_generation(toehold_b, overhang_b, mismatch, type_mismatch):
    mismatch_groups = list(find_range(mismatch))

    base_name = f"{toehold_b}tb_5'{overhang_b}ob"

    full_mismatch_name = str()

    for loc in range(len(mismatch_groups)):
        group = mismatch_groups[loc]
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


def seq_info(fname, trig1="CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"):
    xls = pd.ExcelFile(fname)
    master_list = list()

    for sname in xls.sheet_names:
        dna_seq_df = pd.read_excel(fname, sheet_name=sname)
        loc_df = dna_seq_df.filter(regex="P\d Location", axis=1)
        dna_seq = dna_seq_df["Seqences (5' to 3')"]
        for count in range(len(dna_seq)):
            location = loc_df.iloc[count].dropna()
            for well in location:
                plate = re.findall("P\d", location.index[0])[0]
                plate_loc = "_".join([plate, well])
            seq = dna_seq.iloc[count]
            info = trig_aligment(trig1, seq)
            name = name_generation(*info[:2], *info[3:])
            try:
                dna_trig = trigger(seq, name, *info, plate_loc)
            except UnboundLocalError:
                dna_trig = trigger(seq, name, *info)
            master_list.append(dna_trig)

    trig_df = pd.DataFrame([vars(trig) for trig in master_list])

    return trig_df


if __name__ == "__main__":
    fname = "./dna_sdr/DNA_Strands.xlsx"
    pickle_name = "./dna_sdr/pickles/trig_info.pkl"
    df = seq_info(fname)
    if os.path.exists(pickle_name):
        trig_df = pd.read_pickle(pickle_name)
        if df.equals(trig_df):
            print("Same dataframe - no new file")
        else:
            df.to_pickle(pickle_name)
    else:
        df.to_pickle(pickle_name)
