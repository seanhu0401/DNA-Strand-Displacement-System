"""
Module docstring
"""
import matplotlib.pyplot as plt

from dna_sdr.sequence import seq_utilits
from dna_sdr.sequence.dna_utilits import DNA

marker_dict = {1: ".", 2: "x", 3: "s", 4: "|"}


def mismatched_maker_gen(comparision_dict: dict[int, int]) -> dict:
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
    mismatched_loc = {k: v for k, v in comparision_dict.items() if v != 0}
    mismatched_maker_dict = {}

    for k, v in mismatched_loc.items():
        mismatched_maker_dict[k] = marker_dict[v]

    return mismatched_maker_dict


def seq_comp_plot(seq_1: str, seq_2: str, ax=None):
    """
    xxx

    Parameters
    ----------

    Returns
    -------

    """
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
    T1 = "CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"
    T2 = "CA TAACA CA TCT TT CAATC CA TCT CA CAGCC CA"
    plot = seq_comp_plot(T1, T2)
    plt.show()
