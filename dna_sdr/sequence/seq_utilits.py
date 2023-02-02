from dna_sdr.sequence.dna_utilits import DNA, trigger
import pandas as pd
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


def strand_alignment(
    base_strand: str, incumbent: str, trigger: str, tb_mismatch=False, tb=7
):
    try:
        base_strand = DNA(base_strand)
        incumb = DNA(incumbent)
        trig = DNA(trigger)
    except ValueError:
        raise TypeError("Not all sequence given is DNA sequences")

    bstrand = base_strand.seq
    b_base = base_strand.base_count

    b_incumb = incumb.base_count

    reverse_trig = trig.reverse_complement()
    b_trig = trig.base_count

    counter = mismatch = 0
    base_count = b_incumb + tb
    mismatch_loc = list()

    if tb_mismatch:
        for base in reverse_trig:
            base_at_loc = bstrand[counter]

            if base == base_at_loc and counter + 1 <= b_base:
                pass

            elif base != base_at_loc and counter + 1 <= base_count:
                base_loc_5 = b_trig - counter
                mismatch_loc.append(base_loc_5)
                mismatch += 1

            counter += 1

        ob = b_trig - base_count
        mismatch_loc.sort()

    else:
        last_match_loc = 0
        tb_ob = b_trig - b_incumb  # number of bases for toehold + overhang

        for base in reverse_trig:
            base_at_loc = bstrand[counter]

            if base == base_at_loc and counter <= b_base:
                pass

            elif base != base_at_loc and b_trig - counter > tb_ob:
                base_loc_5 = b_trig - counter
                mismatch_loc.append(base_loc_5)
                mismatch += 1

            elif (
                base != base_at_loc
                and b_trig - counter <= tb_ob
                and counter > last_match_loc
                and last_match_loc < b_trig - counter
            ):
                last_match_loc = counter

            counter += 1

        if last_match_loc != 0:
            toehold = last_match_loc - b_incumb
            ob = tb_ob - toehold

        else:
            tb = tb_ob
            ob = 0

    mismatch_loc.sort()

    return tb, ob, mismatch, mismatch_loc


def name_generation(mismatch, toehold_b, overhang_b):
    mismatch_groups = list(find_range(mismatch))

    base_name = f"{toehold_b}tb_5'{overhang_b}ob"

    for group in mismatch_groups:
        try:
            base_number = group[-1] - group[0] + 1
            mismatch_name = f"_{group[0]}_{base_number}mb"
        except TypeError:
            base_number = 1
            mismatch_name = f"_{group}_{base_number}mb"

    try:
        name = base_name + mismatch_name
    except UnboundLocalError:
        name = base_name

    return name


def seq_info(
    fname,
    incumb="CA CAATC CA TCT CA CCACC CA",
    base_seq="TG GGTGG TG AGA TG GATTG TG AGA TG TG AGA CAT ACA GCG CCG ACC GTA",
):
    xls = pd.ExcelFile(fname)
    master_list = list()

    for sname in xls.sheet_names:
        dna_seq = pd.read_excel(fname, sheet_name=sname)
        print(sname)

        if sname == "no_mtb":
            toehold_bases = 7

        elif sname == "mtb":
            toehold_bases = int(input("toehold base amount:") or "7")

        for seq in dna_seq["Seqences (5' to 3')"]:
            tb, ob, mismatch, mismatch_loc = strand_alignment(
                base_seq, incumb, seq, tb=toehold_bases
            )
            name = name_generation(mismatch_loc, tb, ob)
            dna_trig = trigger(seq, name, tb, ob, mismatch, mismatch_loc)
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
