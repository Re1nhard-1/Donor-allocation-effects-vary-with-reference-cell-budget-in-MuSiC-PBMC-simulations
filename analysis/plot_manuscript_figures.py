from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
OUTPUT = PROJECT.parent / "figures"
ORANGE = "#D55E00"
BLUE = "#0072B2"
GREY = "#B6BCC2"
BUDGETS = [60, 120, 300]

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.6, "lines.linewidth": 0.8,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": "white",
})


def save(fig, name):
    OUTPUT.mkdir(exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUTPUT / f"{name}.{extension}", dpi=450)
    plt.close(fig)
    print(OUTPUT / f"{name}.pdf")


def panel(ax, letter, title):
    ax.set_title(title, loc="left", pad=10)
    ax.text(-0.13, 1.08, letter, transform=ax.transAxes,
            fontweight="bold", fontsize=10, va="bottom")


def budget_contrasts():
    data = pd.read_csv(RESULTS / "music_external_figures/20260917T185640817512Z/external_budget_penalty.source_data.csv")
    fig = plt.figure(figsize=(170 / 25.4, 135 / 25.4))
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.45],
                           left=.095, right=.98, bottom=.17, top=.59, wspace=.75)
    axes = [fig.add_subplot(grid[i]) for i in range(3)]
    for ax, code, title, color, marker in zip(axes[:2], "ab", ["Weighted MuSiC", "Ordinary NNLS"], [ORANGE, BLUE], ["s", "o"]):
        rows = data[data.panel.eq(code)]
        for _, donor in rows[rows.summary.eq("donor")].groupby("held_out", sort=True):
            ax.plot(range(3), donor.set_index("budget").loc[BUDGETS, "value_pp"], color=GREY, linewidth=.65, alpha=.85)
        mean = rows[rows.summary.eq("mean")].set_index("budget").loc[BUDGETS, "value_pp"]
        ax.plot(range(3), mean, color=color, marker=marker, markersize=4, linewidth=1.6, zorder=3)
        ax.axhline(0, color="#383838", linestyle=(0, (3, 3)), linewidth=.7)
        ax.set(xticks=range(3), xticklabels=BUDGETS, ylim=(-.85, 1.1), xlim=(-.16, 2.16),
               yticks=[-.5, 0, .5, 1], xlabel="Cells per type")
        panel(ax, {"a": "b", "b": "c"}[code], title)
    axes[0].set_ylabel("MAE difference (percentage points)")
    ax = axes[2]
    rows = data[data.panel.eq("c") & data.summary.eq("donor")].sort_values("held_out")
    y = np.arange(len(rows))
    ax.hlines(y, 0, rows.value_pp, color=GREY, linewidth=.8)
    ax.scatter(rows.value_pp, y, s=17, marker="s", color=ORANGE, zorder=3)
    mean = data.loc[data.panel.eq("c") & data.summary.eq("mean"), "value_pp"].iloc[0]
    ax.axvline(mean, color=ORANGE, linewidth=1, linestyle=(0, (4, 3)))
    ax.axvline(0, color="#383838", linewidth=.7, linestyle=(0, (3, 3)))
    ax.set(yticks=y, yticklabels=rows.held_out, ylim=(13.65, -.65),
           xlim=(-1.3, .35), xticks=[-1, -.5, 0], xlabel="I(300) − I(60)\n(percentage points)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    panel(ax, "d", "Weighted budget contrast")
    schematic = fig.add_axes([.095, .72, .87, .24])
    schematic.set(xlim=(0, 1), ylim=(0, 1))
    schematic.axis("off")
    schematic.text(-.07, 1, "a", fontweight="bold", fontsize=10, va="top")
    schematic.text(0, 1, "Same three donors; same total cells per type", fontsize=10, va="top")
    colors = ["#0072B2", "#E69F00", "#009E73"]
    for y, counts, label in [(0.58, [20, 20, 20], "Balanced"), (.25, [50, 5, 5], "Unequal")]:
        left = .17
        schematic.text(0, y, label, va="center", fontsize=8)
        for count, color in zip(counts, colors):
            width = .39 * count / 60
            schematic.barh(y, width, left=left, height=.16, color=color, edgecolor="white")
            schematic.text(left + width / 2, y, str(count), ha="center", va="center", color="white", fontsize=7)
            left += width
    schematic.text(.62, .53, "B = 60 → 120 → 300", fontsize=9)
    schematic.text(.62, .22, "Rotate the dominant donor;\nreuse the same target mixtures", fontsize=8)
    schematic.text(.17, .02, "Example at B = 60; colours identify donors", fontsize=7)
    fig.text(.095, .69, "Allocation contrast: 10:1:1 − balanced", fontsize=9, va="top")
    fig.legend(handles=[Line2D([0], [0], color=GREY, label="Held-out donor"),
                        Line2D([0], [0], color="#333333", linewidth=1.6, label="Equal-donor mean (b, c)")],
               loc="lower left", bbox_to_anchor=(.08, .018), frameon=False, ncol=2, fontsize=7.5)
    save(fig, "Figure_1")


def arrangements(method, name, title, color, marker):
    base = RESULTS / "music_arrangement_distribution_v1/20260918T153955745453Z"
    data = pd.read_csv(base / "figure_s4_arrangement_distribution.source_data.csv")
    summary = pd.read_csv(base / "table_s13_distribution_summary.csv")
    data = data[data.method.eq(method)]
    summary = summary[summary.method.eq(method) & summary.source.eq("original_grouping")].set_index("budget")
    donors = sorted(data.held_out.unique())
    fig, axes = plt.subplots(1, 3, figsize=(170 / 25.4, 113 / 25.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.155, right=.985, top=.79, bottom=.24, wspace=.16)
    for ax, budget, letter in zip(axes, BUDGETS, "abc"):
        rows = data[data.budget.eq(budget)]
        for i in range(0, len(donors), 2):
            ax.axhspan(i - .5, i + .5, color="#F2F3F4", zorder=0, linewidth=0)
        ax.axvline(0, color="#383838", linewidth=.8, linestyle=(0, (3, 3)), zorder=1)
        ax.scatter(rows.mae_penalty_pp, rows.plot_y, s=10, marker=marker, c=color,
                   edgecolor="white", linewidth=.25, zorder=2, alpha=.9)
        ax.set(xlim=(-2.5, 3), xticks=[-2, -1, 0, 1, 2, 3], ylim=(13.65, -.65), yticks=range(14), yticklabels=donors)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        panel(ax, letter, f"{budget} cells per type")
        s = summary.loc[budget]
        ax.text(.5, -.15, f"Mean I: {s.signed_mean_pp:+.3f}\nMean |I|: {s.mean_absolute_pp:.3f}",
                ha="center", va="top", transform=ax.transAxes, fontsize=7.5, linespacing=1.5)
    fig.text(.155, .955, title, fontsize=10, va="top")
    fig.text(.155, .90, "12 reference arrangements per held-out donor", fontsize=8, va="top")
    fig.text(.57, .035, "MAE difference: 10:1:1 − balanced (percentage points)", ha="center", fontsize=8)
    save(fig, name)


def diagnostics(main=False):
    data = pd.read_csv(RESULTS / "music_mechanism_review/20260917T211141303956Z/mechanism_diagnostics.source_data.csv")
    codes = "ac" if main else "b"
    fig, axes = plt.subplots(1, len(codes), figsize=((170 if main else 100) / 25.4, 85 / 25.4), squeeze=False)
    fig.subplots_adjust(left=.12 if main else .20, right=.96, top=.80, bottom=.28, wspace=.48)
    titles = {"a": "Profile dispersion", "b": "Between-donor variance", "c": "Fitted gene weights"}
    labels = {"a": "Profile variance across draws", "b": r"Difference in summed Σ ($\times 10^{-6}$)", "c": "Total-variation distance"}
    for index, (ax, code) in enumerate(zip(axes.flat, codes)):
        rows = data[data.panel.eq(code)]
        for _, donor in rows[rows.series.eq("individual_donor")].groupby("donor", sort=True):
            ax.plot(range(3), donor.sort_values("x").value_plotted, color=GREY, linewidth=.7)
        mean = rows[rows.series.eq("equal_donor_mean")].sort_values("x")
        ax.plot(range(3), mean.value_plotted, color=ORANGE, marker="s", markersize=4, linewidth=1.6)
        ax.set(xticks=range(3), xticklabels=mean.x.astype(int), xlim=(-.14, 2.14), ylabel=labels[code])
        ax.set_xlabel("Cells per minor donor and cell type" if code == "a" else "Cells per type")
        panel(ax, chr(97 + index) if main else "", titles[code])
        if code == "a": ax.set_yscale("log")
        if code == "b": ax.axhline(0, color="#383838", linewidth=.7, linestyle=(0, (3, 3)))
        if code == "c": ax.set_ylim(0, .55)
    fig.legend(handles=[Line2D([0], [0], color=GREY, label="Donor summary"),
                        Line2D([0], [0], color=ORANGE, marker="s", markersize=4, label="Equal-donor mean")],
               loc="lower center", bbox_to_anchor=(.5, .015), frameon=False, ncol=2, fontsize=8)
    save(fig, "Figure_3" if main else "Figure_S3")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures", nargs="+", choices=["1", "2", "3", "S3", "S4"], default=["1", "2", "3", "S3", "S4"])
    args = parser.parse_args()
    if "1" in args.figures:
        budget_contrasts()
    if "2" in args.figures:
        arrangements("music_weighted", "Figure_2", "Weighted MuSiC: allocation effects within donors", ORANGE, "s")
    if "3" in args.figures:
        diagnostics(main=True)
    if "S3" in args.figures:
        diagnostics()
    if "S4" in args.figures:
        arrangements("music_nnls", "Figure_S4", "Ordinary NNLS: allocation effects within donors", BLUE, "o")
