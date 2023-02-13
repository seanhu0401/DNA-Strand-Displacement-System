import matplotlib.pyplot as plt
import dna_sdr.sequence.seq_utilits as seq_utilits
from dna_sdr.sequence.dna_utilits import DNA

marker_dict = {1: "o", 2: "x", 3: "s", 4: "|"}


def mismatched_maker_gen(comparision_dict: dict):
    mismatched_loc = {k: v for k, v in comparision_dict.items() if v != 0}
    mismatched_maker_dict = dict()

    for k, v in mismatched_loc.items():
        mismatched_maker_dict[k] = marker_dict[v]

    return mismatched_maker_dict


def seq_comp_plot(seq_1, seq_2, ax=None):
    if ax is None:
        ax = plt.gca()

    comp_dict = seq_utilits.sequence_comparison(seq_1, seq_2)
    seq_tuple = (DNA(seq_1), DNA(seq_2))
    base_loc_lst = list(range(seq_tuple[0].base_count))
    y_temp = [1] * len(base_loc_lst)
    ax.plot(base_loc_lst, y_temp)

    mismatched_maker = mismatched_maker_gen(comp_dict)

    marker_lst = [v for _, v in marker_dict.items()]
    for marker in marker_lst:
        mark_dict = {k: v for k, v in mismatched_maker.items() if v == marker}
        if bool(mark_dict):
            loc_lst = [loc for loc in mark_dict if mark_dict[loc] == marker]
            ax.plot(base_loc_lst, y_temp, marker, markevery=loc_lst, markersize=12)

    ax.set_xlabel("Base number")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax


if __name__ == "__main__":
    t1 = "CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"
    t2 = "CA TAACA CA TCT TT CAATC CA TCT CA CAGCC CA"
    ax = seq_comp_plot(t1, t2)
    plt.show()
