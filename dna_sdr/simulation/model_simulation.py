import pandas as pd
import dna_sdr.simulation.sim_model as model
import dna_sdr.curve_fitting.curve_fitting as cf
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import dna_sdr.visualization.sim_vis as vis

end_plot = True


def model_fitting(pt, parameter, type):
    func_list = list()
    params_list = list()

    for func in dir(model):
        if func.endswith("function"):
            func_list.append(func)
            popt, r_square = cf.general_cf_process(
                getattr(model, func), pt[parameter], pt["mean_{}".format(type)]
            )
            params_list.append((popt, r_square))

    func_para_list = dict(zip(func_list, params_list))

    return func_para_list


def tb_variation_summary(trig_pickle, pt_pickle):
    trig_info_df = pd.read_pickle(trig_pickle)
    pt_out = pd.read_pickle(pt_pickle)
    toehold_trig = trig_info_df[
        (trig_info_df["mismatch"] == 0) & (trig_info_df["overhang"] == 0)
    ].sort_values(by=["toehold"])

    pt_out.columns = ["_".join(col) for col in pt_out.columns.values]
    pt_out.index = ["_".join(col) for col in pt_out.index.values]
    pt_out.reset_index(inplace=True)
    pt_out.rename({"index": "plate_loc"}, axis="columns", inplace=True)

    output = toehold_trig.merge(pt_out, how="inner", on="plate_loc")
    output.drop(
        ["overhang", "mismatch", "mismatch_loc", "seq", "base_count"],
        axis=1,
        inplace=True,
    )

    return output


def mis_loc_variation_summary(trig_pickle, pt_pickle, free_energy_pickle):
    trig_info_df = pd.read_pickle(trig_pickle)
    pt_out = pd.read_pickle(pt_pickle)
    free_energy_df = pd.read_pickle(free_energy_pickle)

    free_energy_df.rename(columns={"Seqences (5' to 3')": "seq"}, inplace=True)
    free_energy_df["seq"] = free_energy_df["seq"].str.lower()
    base = free_energy_df["dG (kcal/mol)"][0]
    free_energy_df["dG (kcal/mol)"] = free_energy_df["dG (kcal/mol)"] / base

    mis_loc = trig_info_df[
        (trig_info_df["mismatch"] == 1) & (trig_info_df["overhang"] == 5)
    ].sort_values(by=["mismatch_loc"])

    pt_out.columns = ["_".join(col) for col in pt_out.columns.values]
    pt_out.index = ["_".join(col) for col in pt_out.index.values]
    pt_out.reset_index(inplace=True)
    pt_out.rename({"index": "plate_loc"}, axis="columns", inplace=True)

    trig_energy_df = mis_loc.merge(free_energy_df, how="inner", on="seq")
    output = trig_energy_df.merge(pt_out, how="inner", on="plate_loc")
    output.drop(
        ["toehold", "overhang", "mismatch", "seq", "base_count"], axis=1, inplace=True
    )

    mis_loc = output["mismatch_loc"]
    int_list = [int(i[0]) for i in mis_loc]
    output["mismatch_loc"] = int_list
    output = output[(output["mismatch_loc"] > 12)]

    return output


if __name__ == "__main__":
    free_energy_pickle = "./dna_sdr/pickles/free_energy.pkl"
    trig_info_pickle = "./dna_sdr/pickles/trig_info.pkl"
    pt = "./dna_sdr/pickles/trig_pt.pkl"
    out = tb_variation_summary(trig_info_pickle, pt)

    result_type = ["plateau", "rate"]
    output = mis_loc_variation_summary(trig_info_pickle, pt, free_energy_pickle)

    # plateau_param = model_fitting(out, "toehold", result_type[0])
    # rate_param = model_fitting(out, "toehold", result_type[1])

    plateau_param = model_fitting(output, "mismatch_loc", result_type[0])
    rate_param = model_fitting(output, "mismatch_loc", result_type[1])

    if end_plot:
        fig = plt.figure(constrained_layout=True)
        gs = GridSpec(2, 1, figure=fig)
        x_inches, y_inches = vis.pt_to_inch(661, 397)

        # ax1 = fig.add_subplot(gs[0])
        # vis.scatter_plot(out, result_type[0], ax=ax1)
        # for func in dir(model):
        #     if func.endswith("function"):
        #         vis.regression_plot(
        #             out, getattr(model, func), plateau_param[func][0], func
        #         )

        # ax1.legend()
        # ax1.set_ylabel("plateau")

        # ax2 = fig.add_subplot(gs[1])
        # vis.scatter_plot(out, result_type[1], ax=ax2)
        # for func in dir(model):
        #     if func.endswith("function"):
        #         vis.regression_plot(
        #             out, getattr(model, func), rate_param[func][0], func
        #         )

        # ax2.set_xlabel("toehold length (bases)")
        # ax2.set_ylabel("rate (1/s)")

        ax1 = fig.add_subplot(gs[0])
        vis.scatter_plot(output, result_type[0], "experimental")
        for func in dir(model):
            if func.endswith("function"):
                vis.regression_plot(
                    output, getattr(model, func), plateau_param[func][0], func
                )

        ax1.legend()
        ax1.set_ylabel("plateau")
        ax1.set_xticks(np.arange(12.0, 34.0, 1.0))

        ax3 = ax1.twinx()
        ax3.scatter(output["mismatch_loc"], output["dG (kcal/mol)"], color="red")
        ax3.set_ylabel("dG ratio", color="red", fontsize=14)
        ax3.tick_params(axis="y", labelcolor="red")

        ax2 = fig.add_subplot(gs[1])
        vis.scatter_plot(output, result_type[1])
        for func in dir(model):
            if func.endswith("function"):
                vis.regression_plot(
                    output, getattr(model, func), rate_param[func][0], func
                )

        ax2.set_xlabel("mismatch location (bases)")
        ax2.set_ylabel("rate (1/s)")
        ax2.set_xticks(np.arange(12.0, 34.0, 1.0))

        ax4 = ax2.twinx()
        ax4.scatter(output["mismatch_loc"], output["dG (kcal/mol)"], color="red")
        ax4.set_ylabel("dG ratio", color="red", fontsize=14)
        ax4.tick_params(axis="y", labelcolor="red")

        plt.show()
