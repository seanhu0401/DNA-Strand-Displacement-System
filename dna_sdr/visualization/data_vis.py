import matplotlib.pyplot as plt
import numpy as np

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


def pt_to_inch(x_pt, y_pt):
    x_inches = x_pt / 72
    y_inches = y_pt / 72
    return x_inches, y_inches


def release_plot_with_error_bar(time, mean, std, ax=None):
    if ax is None:
        ax = plt.gca()

    for y, c in zip(mean, CB_color_cycle):
        ax.scatter(time, mean[y], marker=".", color=c, label=y)
    for y, yerr, c in zip(mean, std, CB_color_cycle):
        ax.errorbar(
            time, mean[y], std[yerr], ls="none", elinewidth=1, capsize=3, color=c
        )

    ax.set_xlabel("Time (mins)", fontsize=12)
    ax.set_ylabel("Normalized % Release", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # legend = ax.legend(loc='lower left', bbox_to_anchor=(1.0, -0.1))
    return ax


def release_plot_with_fit(time, exp_mean, std, fitted_value, ax=None):
    if ax is None:
        ax = plt.gca()

    release_plot_with_error_bar(time, exp_mean, std)

    for y, c in zip(fitted_value, CB_color_cycle):
        ax.plot(time, fitted_value[y], color=c)

    return ax


def final_release_plot(exp_mean, std, ax=None):
    if ax is None:
        ax = plt.gca()

    final_release = exp_mean.iloc[-1, :]
    final_release_std = std.iloc[-1, :]
    label_list = list()
    col_name_list = list(exp_mean.columns)
    for col_name in col_name_list:
        split_name = col_name.split("_")
        label = f"{split_name[0]} ({split_name[1]})"
        label_list.append(label)
    x_pos = np.arange(len(label_list))
    ax.bar(
        x_pos,
        final_release,
        yerr=final_release_std,
        align="center",
        capsize=3,
        color=CB_color_cycle,
    )
    ax.set_xlabel("Type", fontsize=12)
    ax.set_ylabel("Normalized % Release", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(label_list)

    return ax


# TODO: fix the first row of the table to say T1 in the label for Conc. test series
def parameter_table(variable_df, ax=None):
    if ax is None:
        ax = plt.gca()

    ax.axis("tight")
    ax.axis("off")
    table = ax.table(
        cellText=variable_df.values,
        colLabels=variable_df.columns,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(12)

    return table


def params_bar_plot_processing(params_pivot_table):
    cond_list = params_pivot_table.index.values.tolist()
    cond_list.sort(reverse=True)
    cond_list = cond_list[-1:] + cond_list[:-1]

    trig_list = params_pivot_table.columns.values.tolist()
    trig_list = trig_list[-1:] + trig_list[:-1]

    params_pivot_table = params_pivot_table[trig_list]
    params_pivot_table = params_pivot_table.reindex(cond_list)
    return params_pivot_table


def platuea_params_bar_plot(params_pivot_table, ax=None):
    if ax is None:
        ax = plt.gca()

    params_mean = params_pivot_table["mean"]["plateau"]
    params_std = params_pivot_table["std"]["plateau"]

    params_mean = params_bar_plot_processing(params_mean)

    ax = params_mean.T.plot.bar(
        yerr=params_std.T, color=CB_color_cycle, capsize=3, ax=ax
    )

    return ax


def rate_params_bar_plot(params_pivot_table, ax=None):
    if ax is None:
        ax = plt.gca()

    params_mean = params_pivot_table["mean"]["rate_constant"]
    params_std = params_pivot_table["std"]["rate_constant"]

    params_mean = params_bar_plot_processing(params_mean)

    ax = params_mean.T.plot.bar(
        yerr=params_std.T, color=CB_color_cycle, capsize=3, ax=ax
    )

    return ax
