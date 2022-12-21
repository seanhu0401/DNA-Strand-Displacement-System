import regex as re
from dna_class import DNA

base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}
comp_dict = {"t": 1, "g": 2, "c": 3, "a": 4}


def remove(seq: str):
    return seq.replace(" ", "")


def is_dna(seq: str):
    lower_case_seq = remove(seq).lower()
    return bool(re.match("^[atgc]+$", lower_case_seq))


def reverse_seq(seq: str):
    return seq[::-1]


# *Think of better name for the method
# This function is currently used to compare to sequence, mainly the triggers, to identify the location(s) that has/have mismatch(es)
# They should be equal length in bases if there is no modification to the toehold + dangles - not applicable with those modification currently
def sequence_comparison(seq_1: DNA, seq_2: DNA):
    comparison_lst = []
    base_loc_lst = []
    counter = 0

    if seq_1.base_count != seq_2.base_count:
        raise Exception("The two sequences are not equal length")

    while counter <= seq_1.base_count:
        comp_val = 0
        if seq_1.seq[counter] == seq_2.seq[counter]:
            comparison_lst.append(comp_val)
        elif seq_1.seq[counter] != seq_2.seq[counter]:
            comp_val = comp_dict[seq_2[counter]]
            comparison_lst.append(comp_val)
        base_loc_lst.append(counter + 1)
        counter += 1

    seq_comp_dict = dict(zip(base_loc_lst, comparison_lst))
    return seq_comp_dict
