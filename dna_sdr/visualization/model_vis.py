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

CB_color_cycle = [
    "#377eb8",
    "#ff7f00",
    "#4daf4a",
    "#f781bf",
    "#a65628",
    "#984ea3",
    "#999999",
    "#e41a1c",
    "#dede00",
]


def fit_plot(test, fit_cond, x, result, condition, color="#377eb8", ax=None):
    if ax is None:
        ax = plt.gca()

    data = result.data
    uncertainty = result.eval_uncertainty(sigma=3)
    error = 1 / result.weights

    if test == "Conc":
        ax.errorbar(
            x,
            data,
            error,
            elinewidth=1,
            capsize=3,
            fmt="o",
            label="{}%".format(condition.split("_")[-1]),
            color=color,
        )
        ax.plot(x, result.best_fit, "--", color="#ff7f00")
        ax.fill_between(
            x,
            result.best_fit - uncertainty,
            result.best_fit + uncertainty,
            color="#ABABAB",
        )
    elif test == "screen":
        ax.errorbar(
            x, data, error, elinewidth=1, capsize=3, fmt="o", label="exp.", color=color
        )
        ax.fill_between(
            x,
            result.best_fit - uncertainty,
            result.best_fit + uncertainty,
            color="#ABABAB",
            label="3-$\sigma$ uncertainty band",
        )
        ax.plot(x, result.best_fit, "--", label=fit_cond, color="#ff7f00")
        ax.fill_between(
            x,
            result.best_fit - uncertainty,
            result.best_fit + uncertainty,
            color="#ABABAB",
            label="3-$\sigma$ uncertainty band",
        )
        ax.text(
            0.95,
            0.1,
            "r$^2$ = {:2f}".format(result.rsquared),
            verticalalignment="bottom",
            horizontalalignment="right",
            transform=ax.transAxes,
            fontsize=15,
        )
    ax.set_ylabel("Release (nM)")
    ax.set_xlabel("Time (min)")
    return ax


if __name__ == "__main__":
    os.chdir("./dna_sdr/pickles/")
    time = time_list_generation(60)
    test = "Conc"

    with open("{}_one_phase_result.pkl".format(test), "rb") as one_phase_assoc:
        one_phase = pickle.load(one_phase_assoc)

    with open("{}_kinetic_result.pkl".format(test), "rb") as kin:
        kine = pickle.load(kin)

    if test == "Conc":
        filter_one_phase_lst = list()
        filter_kine_lst = list()
        x = ["_".join(condition.split("_")[:2]) for condition in one_phase]
        con = [*set(x)]
        for item in con:
            filtered_one_phase = {
                k: v for k, v in one_phase.items() if (item in k and k != "P0_T1")
            }
            filter_one_phase_lst.append(filtered_one_phase)
            filtered_kin = {
                k: v for k, v in kine.items() if (item in k and k != "P0_T1")
            }
            filter_kine_lst.append(filtered_kin)

        for set in range(len(filter_one_phase_lst)):
            fig, (ax1, ax2) = plt.subplots(1, 2)
            counter = 0
            for condition in filter_one_phase_lst[set].keys():
                one_phase_result = filter_one_phase_lst[set][condition]
                kine_result = filter_kine_lst[set][condition]
                fit_plot(
                    test,
                    "one phase assoc.",
                    time,
                    one_phase_result,
                    condition,
                    color=CB_color_cycle[counter],
                    ax=ax1,
                )
                fit_plot(
                    test,
                    "kinetic",
                    time,
                    kine_result,
                    condition,
                    color=CB_color_cycle[counter],
                    ax=ax2,
                )
                ax1.set_title("One phase association")
                ax2.set_title("kinetic")
                ax1.legend()
                ax2.legend()
                counter += 1

            fig.suptitle("_".join(condition.split("_")[0:2]))
            fig.set_figheight(9)
            fig.set_figwidth(16)
            # fig.savefig("./{}/{}.svg".format(test, condition), format="svg")
            plt.show()
            # plt.close()

    elif test == "screen":
        for condition in one_phase:
            result_one_phase = one_phase[condition]

            result_kin = kine[condition]

            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1 = fit_plot(test, "one phase assoc.", time, result_one_phase, ax=ax1)
            ax2 = fit_plot(test, "kinetic", time, result_kin, ax=ax2)

            fig.suptitle(condition)
            fig.set_figheight(9)
            fig.set_figwidth(16)
            ax1.legend()
            ax2.legend()
            # fig.savefig("./{}/{}.svg".format(test, condition), format="svg")
            # plt.show()
            # plt.close()
