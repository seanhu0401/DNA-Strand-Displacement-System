import matplotlib.pyplot as plt
import dna_sdr.sequence.seq_utilits as seq_utilits
from dna_sdr.sequence.dna_utilits import DNA


def seq_comp_plot(seq_1, seq_2, ax=None):
    if ax is None:
        ax = plt.gca()

    seq_tuple = (DNA(seq_1), DNA(seq_2))
    comp_dict = seq_utilits.sequence_comparison(*seq_tuple)
    base_loc_lst = list(range(seq_tuple[0].base_count))

    ax.set_xlabel("Base number")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


if __name__ == "__main__":
    seq1 = "TCT CA CAATC CA TCT"
    seq2 = "TCT CA TAATC CA TCT"
    seq_comp_plot(seq1, seq2)
