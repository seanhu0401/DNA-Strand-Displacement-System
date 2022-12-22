import matplotlib.pyplot as plt
import dna_sdr.sequence.seq_utilits as seq_utilits
from dna_sdr.sequence.dna_utilits import DNA

maker_dict = {1: "o", 2: "x", 3: "s", 4: "*"}


def mismatched_maker_gen(comparision_dict: dict):
    mismatched_loc = {k: v for k, v in comparision_dict.items() if v != 0}
    mismatched_loc_lst = list()
    marker_lst = list()

    for k, v in mismatched_loc.items():
        mismatched_loc_lst.append(k)
        marker = maker_dict[v]
        marker_lst.append(marker)

    return mismatched_loc_lst, marker_lst


def seq_comp_plot(seq_1, seq_2, ax=None):
    if ax is None:
        ax = plt.gca()

    seq_tuple = (DNA(seq_1), DNA(seq_2))
    comp_dict = seq_utilits.sequence_comparison(*seq_tuple)
    base_loc_lst = list(range(seq_tuple[0].base_count))
    y_temp = [1] * len(base_loc_lst)
    ax.plot(base_loc_lst, y_temp)

    marker, mismatch_loc = mismatched_maker_gen(comp_dict)

    ax.set_xlabel("Base number")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax
