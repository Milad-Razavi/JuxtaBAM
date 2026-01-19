#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

#############################################
# Global style defaults
#############################################

MAX_FIG_WIDTH = 40
MAX_FIG_HEIGHT = 30
SAVE_DPI = 350

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 20,
})

color_cmap = LinearSegmentedColormap.from_list(
    "lightgray_red",
    ["#f5f5f5", "#fee5d9", "#fcae91", "#fb6a4a", "#8B0000"]
).reversed()

#############################################
def scaled_fonts(n_labels):
    base = 18
    scale = min(1.0, 60 / max(1, n_labels))

    tick_fs = max(10, int(base * scale))
    label_fs = max(12, int(22 * scale))
    legend_fs = max(10, int(18 * scale))
    title_fs = max(14, int(26 * scale))
    cbar_fs = max(10, int(18 * scale))
    return tick_fs, label_fs, legend_fs, title_fs, cbar_fs


def compute_bubble_sizes(vals, size_min=20, size_max=350):
    vmin = vals.min()
    vmax = vals.max()

    if vmax > vmin:
        vals_norm = (vals - vmin) / (vmax - vmin)
    else:
        vals_norm = np.zeros_like(vals)

    sizes = size_min + (1.0 - vals_norm) * (size_max - size_min)
    return sizes, vmin, vmax


def plot_heatmap(df, a1, a2, outdir, min_annot_percent):

    sub = df[(df["assay1"] == a1) & (df["assay2"] == a2)]
    if sub.empty:
        return

    heat = sub.pivot(index="sample1", columns="sample2",
                    values="percent_mismatch")

    annot = heat.copy().astype(str)

    for i in annot.index:
        for j in annot.columns:
            row = sub[(sub["sample1"] == i) &
                    (sub["sample2"] == j)]
            if row.empty:
                annot.at[i, j] = ""
                continue

            v = row["percent_mismatch"].iloc[0]
            b1 = row["biosample1"].iloc[0]
            b2 = row["biosample2"].iloc[0]

            if i == j or b1 == b2:
                annot.at[i, j] = f"{v:.1f}"
            elif pd.notna(v) and v < min_annot_percent:
                annot.at[i, j] = f"{v:.1f}"
            else:
                annot.at[i, j] = ""

    nx = heat.shape[1]
    ny = heat.shape[0]
    n = max(nx, ny)

    width = min(MAX_FIG_WIDTH, max(12, 0.45 * nx))
    height = min(MAX_FIG_HEIGHT, max(10, 0.45 * ny))
    plt.figure(figsize=(width, height))

    tick_fs, label_fs, legend_fs, title_fs, cbar_fs = scaled_fonts(n)

    sns.heatmap(
        heat,
        annot=annot,
        fmt="",
        cmap=color_cmap,
        linewidths=0.5,
        linecolor="gray",
        cbar_kws={"label": "Genotype Mismatch (%)"}
    )

    cbar = plt.gca().collections[0].colorbar
    cbar.ax.tick_params(labelsize=cbar_fs)
    cbar.set_label("Genotype Mismatch (%)", size=label_fs)

    plt.xticks(rotation=90, fontsize=tick_fs)
    plt.yticks(rotation=0, fontsize=tick_fs)

    plt.title(f"{a1} vs {a2} - SNV Mismatch Rate", fontsize=title_fs)

    plt.xlabel(f"{a2} samples", fontsize=label_fs)
    plt.ylabel(f"{a1} samples", fontsize=label_fs)

    plt.tight_layout()
    plt.savefig(outdir / f"heatmap_{a1}_vs_{a2}.png", dpi=SAVE_DPI)
    plt.savefig(outdir / f"heatmap_{a1}_vs_{a2}.pdf")
    plt.close()


def plot_bubble_heatmap(df, a1, a2, outdir):
    sub = df[(df["assay1"] == a1) & (df["assay2"] == a2)]
    if sub.empty:
        return

    heat = sub.pivot(index="sample1", columns="sample2",
                    values="percent_mismatch")

    samples_y = heat.index.tolist()
    samples_x = heat.columns.tolist()

    nx = len(samples_x)
    ny = len(samples_y)
    n = max(nx, ny)

    width = min(MAX_FIG_WIDTH, max(12, 0.45 * nx))
    height = min(MAX_FIG_HEIGHT, max(10, 0.45 * ny))
    plt.figure(figsize=(width, height))

    tick_fs, label_fs, legend_fs, title_fs, cbar_fs = scaled_fonts(n)

    xs, ys, vals = [], [], []
    for iy, y in enumerate(samples_y):
        for ix, x in enumerate(samples_x):
            v = heat.at[y, x]
            if pd.notna(v):
                xs.append(ix)
                ys.append(iy)
                vals.append(v)

    vals = np.array(vals)

    DOT_SIZE_MIN = 20
    DOT_SIZE_MAX = 350

    sizes, vmin, vmax = compute_bubble_sizes(vals, DOT_SIZE_MIN, DOT_SIZE_MAX)

    sc = plt.scatter(xs, ys, c=vals, s=sizes,
                    cmap=color_cmap, alpha=0.85, edgecolors="none")

    legend_handles = [
        plt.scatter([], [], s=DOT_SIZE_MAX, color=color_cmap(0), alpha=0.85),
        plt.scatter([], [], s=DOT_SIZE_MIN, color=color_cmap(0), alpha=0.85),
    ]

    plt.gcf().legend(
        legend_handles,
        [f"Min mismatch ({vmin:.1f}%)", f"Max mismatch ({vmax:.1f}%)"],
        loc="lower left",
        bbox_to_anchor=(0.01, 0.01),
        frameon=True,
        fontsize=legend_fs
    )

    plt.xticks(range(nx), samples_x, rotation=90, fontsize=tick_fs)
    plt.yticks(range(ny), samples_y, fontsize=tick_fs)

    plt.xlabel(f"{a2} samples", fontsize=label_fs)
    plt.ylabel(f"{a1} samples", fontsize=label_fs)
    plt.title(f"{a1} vs {a2} - SNV Mismatch Rate", fontsize=title_fs)

    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")

    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=cbar_fs)
    cbar.set_label("Genotype Mismatch (%)", size=label_fs)

    plt.tight_layout()
    plt.savefig(outdir / f"bubble_heatmap_{a1}_vs_{a2}.png", dpi=SAVE_DPI)
    plt.savefig(outdir / f"bubble_heatmap_{a1}_vs_{a2}.pdf")
    plt.close()


def plot_all_vs_all_bubble(df, outdir):

    heat = df.pivot(index="sample1", columns="sample2",
                    values="percent_mismatch")

    samples_y = heat.index.tolist()
    samples_x = heat.columns.tolist()

    nx = len(samples_x)
    ny = len(samples_y)
    n = max(nx, ny)

    width = min(MAX_FIG_WIDTH, max(14, 0.45 * nx))
    height = min(MAX_FIG_HEIGHT, max(12, 0.45 * ny))
    plt.figure(figsize=(width, height))

    tick_fs, label_fs, legend_fs, title_fs, cbar_fs = scaled_fonts(n)

    xs, ys, vals = [], [], []
    for iy, y in enumerate(samples_y):
        for ix, x in enumerate(samples_x):
            v = heat.at[y, x]
            if pd.notna(v):
                xs.append(ix)
                ys.append(iy)
                vals.append(v)

    vals = np.array(vals)

    DOT_SIZE_MIN = 20
    DOT_SIZE_MAX = 350

    sizes, vmin, vmax = compute_bubble_sizes(vals, DOT_SIZE_MIN, DOT_SIZE_MAX)

    sc = plt.scatter(xs, ys, c=vals, s=sizes,
                    cmap=color_cmap, alpha=0.85, edgecolors="none")

    legend_handles = [
        plt.scatter([], [], s=DOT_SIZE_MAX, color=color_cmap(0), alpha=0.85),
        plt.scatter([], [], s=DOT_SIZE_MIN, color=color_cmap(0), alpha=0.85),
    ]

    plt.gcf().legend(
        legend_handles,
        [f"Min mismatch ({vmin:.1f}%)", f"Max mismatch ({vmax:.1f}%)"],
        loc="lower left",
        bbox_to_anchor=(0.01, 0.01),
        frameon=True,
        fontsize=legend_fs
    )

    plt.xticks(range(nx), samples_x, rotation=90, fontsize=tick_fs)
    plt.yticks(range(ny), samples_y, fontsize=tick_fs)

    plt.xlabel("Samples", fontsize=label_fs)
    plt.ylabel("Samples", fontsize=label_fs)
    plt.title("All samples vs All samples - SNV Mismatch Rate",
            fontsize=title_fs)

    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")

    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=cbar_fs)
    cbar.set_label("Genotype Mismatch (%)", size=label_fs)

    plt.tight_layout()
    plt.savefig(outdir / "bubble_heatmap_all_vs_all.png", dpi=SAVE_DPI)
    plt.savefig(outdir / "bubble_heatmap_all_vs_all.pdf")
    plt.close()


def main():
    #############################################
    # Arguments
    #############################################
    parser = argparse.ArgumentParser()

    parser.add_argument("--genotype_mismatch_rates", required=True)
    parser.add_argument("--input_info", required=True)
    parser.add_argument("--outdir", default="output_plots")
    parser.add_argument("--min_annot_percent", type=float, default=10)
    parser.add_argument("--no_line_annotations", action="store_true")

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)
    #############################################
    # Load & merge metadata
    #############################################
    comparisons = pd.read_csv(args.genotype_mismatch_rates, sep="\t")
    manifest = pd.read_csv(args.input_info, sep="\t")

    comparisons = comparisons.merge(
        manifest[["sample_id", "biosample", "assay"]],
        left_on="sample1",
        right_on="sample_id",
        how="left"
    ).rename(columns={
        "biosample": "biosample1",
        "assay": "assay1"
    }).drop(columns=["sample_id"])

    comparisons = comparisons.merge(
        manifest[["sample_id", "biosample", "assay"]],
        left_on="sample2",
        right_on="sample_id",
        how="left"
    ).rename(columns={
        "biosample": "biosample2",
        "assay": "assay2"
    }).drop(columns=["sample_id"])

    comparisons["percent_mismatch"] = pd.to_numeric(
        comparisons["percent_mismatch"], errors="coerce"
    )

    #############################################
    # Make symmetric sample table
    #############################################
    rev = comparisons.copy()
    for base in ["sample", "biosample", "assay"]:
        rev[f"{base}1"], rev[f"{base}2"] = comparisons[f"{base}2"], comparisons[f"{base}1"]

    comparisons = pd.concat([comparisons, rev], ignore_index=True)
    comparisons = comparisons.drop_duplicates(
        subset=["sample1", "sample2"], keep="first"
    ).reset_index(drop=True)

    assays = sorted(comparisons["assay1"].dropna().unique())

    #############################################
    # Number of SNV sites histograms
    #############################################
    for a1, a2 in product(assays, repeat=2):
        
        sub = comparisons[(comparisons["assay1"] == a1) &
                        (comparisons["assay2"] == a2)]
        
        if sub.empty:
            continue

        same = sub[sub["biosample1"] == sub["biosample2"]]

        vals_all = sub["total_sites"].dropna()
        vals_same = same["total_sites"].dropna()

        if len(vals_all) == 0:
            continue

        # dynamic figure width based on sample count
        n = len(vals_all)
        width = min(MAX_FIG_WIDTH, max(10, 0.35 * n))
        plt.figure(figsize=(width, 8))


        tick_fs, label_fs, legend_fs, title_fs, cbar_fs = scaled_fonts(n)

        plt.hist(vals_all, bins="auto",
                alpha=0.45, color="C1", edgecolor="black",
                label="All pairs")

        if len(vals_same) > 0:
            plt.hist(vals_same, bins="auto",
                    alpha=0.45, color="C0", edgecolor="black",
                    label="Same biosample pairs")

        plt.xlabel("Number of genotyped SNV sites", fontsize=label_fs)
        plt.ylabel("Count", fontsize=label_fs)
        plt.title(f"SNV site distribution: {a1} vs {a2}", fontsize=title_fs)

        plt.xticks(fontsize=tick_fs)
        plt.yticks(fontsize=tick_fs)

        plt.legend(fontsize=legend_fs)

        plt.tight_layout()
        plt.savefig(outdir / f"snv_sites_{a1}_vs_{a2}.png", dpi=SAVE_DPI)
        plt.savefig(outdir / f"snv_sites_{a1}_vs_{a2}.pdf")
        plt.close()


    #############################################
    # LINE PLOTS
    #############################################
    for a1, a2 in product(assays, repeat=2):

        if a1 == a2:
            continue

        sub = comparisons[(comparisons["assay1"] == a1) &
                        (comparisons["assay2"] == a2)]

        if sub.empty:
            continue

        same = sub[sub["biosample1"] == sub["biosample2"]]
        diff = sub[sub["biosample1"] != sub["biosample2"]]

        self_min = same.groupby("biosample1")["percent_mismatch"].min()
        nonself_min = diff.groupby("biosample1")["percent_mismatch"].min()

        bios = sorted(set(self_min.index) | set(nonself_min.index))
        if not bios:
            continue

        labels = [f"{a1}_{b}" for b in bios]

        y1 = self_min.reindex(bios)
        y2 = nonself_min.reindex(bios)

        x = range(len(bios))

        # dynamic sizing
        n = len(bios)
        width = min(MAX_FIG_WIDTH, max(12, 0.5 * n))
        plt.figure(figsize=(width, 10))

        tick_fs, label_fs, legend_fs, title_fs, cbar_fs = scaled_fonts(n)

        plt.plot(x, y1.values, marker="o", color="C0",
                label="Same biosample")

        plt.plot(x, y2.values, marker="o", color="C1",
                label="Other biosamples (min)")

        # annotate potential mismatches on the lineplot, if enabled
        if not args.no_line_annotations:
            for idx, b in enumerate(bios):
                if pd.isna(y1[b]) or pd.isna(y2[b]):
                    continue
                if y1[b] > y2[b]:
                    row = diff[(diff["biosample1"] == b)].sort_values(
                        "percent_mismatch"
                    ).head(1)
                    if not row.empty:
                        label = f"{row['assay2'].iloc[0]}_{row['biosample2'].iloc[0]}"
                        plt.text(idx, y2[b], label,
                                fontsize=max(9, int(tick_fs * 0.9)),
                                ha="left", va="bottom", color="C1")

        plt.xticks(x, labels, rotation=90, fontsize=tick_fs)
        plt.yticks(fontsize=tick_fs)

        plt.ylabel("Percent mismatch (%)", fontsize=label_fs)
        plt.xlabel("", fontsize=label_fs)
        plt.title(f"{a1} vs {a2}: SNV Mismatch Rate", fontsize=title_fs)
        plt.legend(fontsize=legend_fs)

        plt.tight_layout()
        plt.savefig(outdir / f"lineplot_{a1}_vs_{a2}.png", dpi=SAVE_DPI)
        plt.savefig(outdir / f"lineplot_{a1}_vs_{a2}.pdf")

        plt.close()


    #############################################
    # ALL-vs-ALL HEATMAP
    #############################################
    heat = comparisons.pivot(
        index="sample1", columns="sample2",
        values="percent_mismatch"
    )

    annot = heat.copy().astype(str)

    for i in annot.index:
        for j in annot.columns:
            row = comparisons[(comparisons["sample1"] == i) &
                            (comparisons["sample2"] == j)]

            if row.empty:
                annot.at[i, j] = ""
                continue

            v = row["percent_mismatch"].iloc[0]
            b1 = row["biosample1"].iloc[0]
            b2 = row["biosample2"].iloc[0]

            if i == j or b1 == b2:
                annot.at[i, j] = f"{v:.1f}"
            elif pd.notna(v) and v < args.min_annot_percent:
                annot.at[i, j] = f"{v:.1f}"
            else:
                annot.at[i, j] = ""

    n = heat.shape[0]
    width = min(MAX_FIG_WIDTH, max(14, 0.45 * n))
    height = min(MAX_FIG_HEIGHT, max(12, 0.45 * n))
    plt.figure(figsize=(width, height))

    tick_fs, label_fs, legend_fs, title_fs, cbar_fs = scaled_fonts(n)

    sns.heatmap(
        heat,
        annot=annot,
        fmt="",
        cmap=color_cmap,
        linewidths=0.3,
        linecolor="gray",
        cbar_kws={"label": "Genotype Mismatch (%)"}
    )

    cbar = plt.gca().collections[0].colorbar
    cbar.ax.tick_params(labelsize=cbar_fs)
    cbar.set_label("Genotype Mismatch (%)", size=label_fs)

    plt.xticks(rotation=90, fontsize=tick_fs)
    plt.yticks(rotation=0, fontsize=tick_fs)

    plt.title("All samples vs All samples - SNV Mismatch Rate",
            fontsize=title_fs)

    plt.xlabel("")
    plt.ylabel("")

    plt.tight_layout()
    plt.savefig(outdir / "heatmap_all_vs_all.png", dpi=SAVE_DPI)
    plt.savefig(outdir / "heatmap_all_vs_all.pdf")
    plt.close()

    #############################################
    # Run plots
    #############################################
    for a1, a2 in product(assays, repeat=2):
        plot_heatmap(comparisons, a1, a2, outdir, args.min_annot_percent)
        plot_bubble_heatmap(comparisons, a1, a2, outdir)

    plot_all_vs_all_bubble(comparisons, outdir)

if __name__ == "__main__":
    main()
