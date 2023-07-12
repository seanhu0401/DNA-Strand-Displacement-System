import os, pickle
import numpy as np
import matplotlib.pyplot as plt
from dna_sdr.fitting.fit import one_phase_association, kin_fit
from dna_sdr.data_process.list_generation import time_list_generation

plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("legend", fontsize=14)
plt.rc("axes", labelsize=14)
plt.rc("axes", titlesize=14)


if __name__ == "__main__":
    os.chdir("./dna_sdr/pickles/")
    time = time_list_generation(60)
    test = "Conc"

    with open("{}_one_phase_result.pkl".format(test), "rb") as one_phase_assoc:
        one_phase = pickle.load(one_phase_assoc)

    with open("{}_kinetic_result.pkl".format(test), "rb") as kin:
        kine = pickle.load(kin)

    for condition in one_phase:
        result_one_phase = one_phase[condition]
        data_one_phase = result_one_phase.data
        dely_one_phase = result_one_phase.eval_uncertainty(sigma=3)

        result_kin = kine[condition]
        data_kin = result_kin.data
        dely_kin = result_kin.eval_uncertainty(sigma=3)

        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.errorbar(
            time,
            data_one_phase,
            1 / result_one_phase.weights,
            elinewidth=1,
            capsize=3,
            fmt="o",
            label="exp.",
        )
        ax1.plot(
            time,
            result_one_phase.best_fit,
            "--",
            label="one phase association",
            color="r",
        )
        ax1.fill_between(
            time,
            result_one_phase.best_fit - dely_one_phase,
            result_one_phase.best_fit + dely_one_phase,
            color="#ABABAB",
            label="3-$\sigma$ uncertainty band",
        )
        ax1.text(
            0.95,
            0.2,
            "r$^2$ = {:2f}".format(result_one_phase.rsquared),
            verticalalignment="bottom",
            horizontalalignment="right",
            transform=ax1.transAxes,
            fontsize=15,
        )
        ax1.set_ylabel("Release (nM)")

        ax2.errorbar(
            time,
            data_kin,
            1 / result_kin.weights,
            elinewidth=1,
            capsize=3,
            fmt="o",
            label="exp.",
        )
        ax2.plot(time, result_kin.best_fit, "--", label="kinetic")
        ax2.fill_between(
            time,
            result_kin.best_fit - dely_kin,
            result_kin.best_fit + dely_kin,
            color="#ABABAB",
            label="3-$\sigma$ uncertainty band",
        )
        ax2.text(
            0.95,
            0.2,
            "r$^2$ = {:2f}".format(result_kin.rsquared),
            verticalalignment="bottom",
            horizontalalignment="right",
            transform=ax2.transAxes,
            fontsize=15,
        )
        ax2.set_ylabel("Release (nM)")

        for ax in fig.get_axes():
            ax.set_xlabel("Time (min)")

        fig.suptitle(condition)
        fig.set_figheight(9)
        fig.set_figwidth(16)
        ax1.legend()
        ax2.legend()
        fig.savefig("./{}/{}.svg".format(test, condition), format="svg")
        plt.close()
