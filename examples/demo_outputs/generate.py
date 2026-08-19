"""Generate independent GenomeCanvas demonstration figures.

The figures use deterministic synthetic data and demonstrate the same broad
publication behaviors as common genomic figure demos: faceting, isoforms,
interval evidence, explicit signal ratios, highlights, and zoom panels. This
script is original example code and does not reproduce another project's
implementation or image templates.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from genomecanvas.core.intervals import GenomicInterval
from genomecanvas.layout.interval_packing import pack_intervals
from genomecanvas.render.plotting import make_grid, plot_signal, save_figure, style_axis

OUT = Path(__file__).parent
X = np.arange(0, 1000)
SIGNAL_A = 0.2 + np.exp(-((X - 270) / 75) ** 2) + 0.7 * np.exp(-((X - 690) / 55) ** 2)
SIGNAL_B = 0.15 + 0.55 * np.exp(-((X - 300) / 90) ** 2) + 1.0 * np.exp(-((X - 700) / 60) ** 2)
EXONS = [(80, 210), (300, 390), (530, 650), (760, 920)]


def draw_gene_model(axis: plt.Axes, *, reverse: bool = False, labels: bool = True) -> None:
    """Draw an independent transcript model beneath a signal panel."""
    axis.plot([0.08, 0.9], [0.62, 0.62], color="#444444", linewidth=0.8)
    if reverse:
        axis.annotate(
            "", xy=(0.08, 0.62), xytext=(0.9, 0.62),
            arrowprops={"arrowstyle": "->", "color": "#777"},
        )
    else:
        axis.annotate(
            "", xy=(0.9, 0.62), xytext=(0.08, 0.62),
            arrowprops={"arrowstyle": "->", "color": "#777"},
        )
    exons = [(0.12, 0.24), (0.38, 0.45), (0.62, 0.72), (0.78, 0.88)]
    for start, end in exons:
        axis.add_patch(plt.Rectangle((start, 0.48), end - start, 0.28, color="#222222"))
    if labels:
        axis.text(0.02, 0.62, "Gene model", fontsize=7, va="center", ha="right")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    axis.set_xticks([])
    axis.spines[["top", "right", "left", "bottom"]].set_visible(False)


def draw_cytoband(axis: plt.Axes, *, marker: float = 0.12) -> None:
    """Draw a schematic chromosome ideogram with a locus marker."""
    bands = [
        (0.0, 0.1, "#111111"), (0.1, 0.2, "#bdbdbd"), (0.2, 0.3, "#f7f7f7"),
        (0.3, 0.4, "#8c8c8c"), (0.4, 0.5, "#f7f7f7"), (0.5, 0.62, "#111111"),
        (0.62, 0.72, "#bdbdbd"), (0.72, 0.82, "#f7f7f7"), (0.82, 1.0, "#111111"),
    ]
    for start, end, color in bands:
        axis.add_patch(
            plt.Rectangle((start, 0.3), end - start, 0.4, color=color, ec="#555555", lw=0.3)
        )
    axis.plot([0, 1], [0.22, 0.22], color="#222222", linewidth=0.7)
    axis.scatter([marker], [0.1], marker="v", s=70, color="#e6614f", zorder=3)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    axis.set_xticks([])
    axis.spines[["top", "right", "left", "bottom"]].set_visible(False)


def plot_faceted() -> None:
    """Four locus facets, each with four aligned sample signal tracks."""
    loci = ["Locus-A", "Locus-B", "Locus-C", "Locus-D"]
    fig, axes = make_grid(6, 4, sharex="col", figsize=(13, 10),
                          height_ratios=[2, 2, 2, 2, 0.8, 0.45], constrained_layout=True)
    for column, locus in enumerate(loci):
        rng = np.random.default_rng(column + 7)
        local_x = np.linspace(0, 1, 180)
        for row, label in enumerate(["Input rep1", "Input rep2", "IP rep1", "IP rep2"]):
            centers = rng.uniform(0.15, 0.9, size=4)
            values = sum(
                np.exp(-((local_x - center) / rng.uniform(0.015, 0.06)) ** 2)
                for center in centers
            )
            values *= 55 + 25 * row
            color = "#d0d0d0" if row < 2 else "#606060"
            axis = axes[row, column]
            plot_signal(axis, local_x, values, color=color, alpha=1.0, linewidth=0)
            axis.set_ylim(0, max(100, float(values.max()) * 1.1))
            axis.set_yticks([0, int(axis.get_ylim()[1] / 2)])
            axis.tick_params(labelsize=7, length=2)
            axis.spines[["top", "right", "left"]].set_color("#d7d7d7")
            if column == 0:
                axis.set_ylabel(label, fontsize=8)
            if row == 0:
                axis.set_title(
                    f"{locus} | chr1:10,000-{11_000 + column * 900:,}",
                    fontsize=9, weight="bold",
                )
        draw_gene_model(axes[4, column], reverse=column % 2 == 1)
        draw_cytoband(axes[5, column], marker=0.1 + column * 0.18)
    fig.suptitle("Synthetic locus comparison · four aligned facets", fontsize=13, weight="bold")
    save_figure(fig, OUT, "p01_faceted")


def plot_isoforms() -> None:
    """Four locus facets with signal rows and individually packed transcripts."""
    fig, axes = make_grid(6, 4, figsize=(13, 11), sharex="col",
                          height_ratios=[2, 2, 2, 2, 1.6, 0.45], constrained_layout=True)
    for column, gene in enumerate(["Gene-A", "Gene-B", "Gene-C", "Gene-D"]):
        local_x = np.linspace(0, 1, 180)
        for row in range(4):
            values = (1 + row * 0.5) * (
                np.exp(-((local_x - 0.28 - row * 0.02) / 0.07) ** 2)
                + 0.6 * np.exp(-((local_x - 0.7) / 0.045) ** 2)
            )
            axis = axes[row, column]
            plot_signal(axis, local_x, values, color="#bdbdbd" if row < 2 else "#4d4d4d",
                        alpha=1.0, linewidth=0)
            axis.set_yticks([])
            axis.spines[["top", "right", "left"]].set_color("#d7d7d7")
            if column == 0:
                axis.set_ylabel(["Input rep1", "Input rep2", "IP rep1", "IP rep2"][row], fontsize=8)
            if row == 0:
                axis.set_title(
                    f"{gene} | chr1:20,000-{21_000 + column * 500:,}",
                    fontsize=9, weight="bold",
                )
        axis = axes[4, column]
        for transcript, (start, end) in enumerate([(0.08, 0.92), (0.18, 0.78), (0.3, 0.96)]):
            y = 0.8 - transcript * 0.28
            axis.plot([start, end], [y, y], color="#555555", linewidth=0.7)
            exons = [(start, start + 0.1), (start + 0.25, start + 0.38), (end - 0.1, end)]
            for exon_start, exon_end in exons:
                axis.add_patch(
                    plt.Rectangle(
                        (exon_start, y - 0.07), exon_end - exon_start, 0.14, color="#222222"
                    )
                )
            axis.text(0.02, y, f"ENST{column + 1:02d}{transcript + 1:02d}", fontsize=7, va="center")
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_yticks([])
        axis.set_xlabel("coordinate", fontsize=8)
        axis.spines[["top", "right", "left"]].set_visible(False)
        draw_cytoband(axes[5, column], marker=0.16 + column * 0.17)
    fig.suptitle("Synthetic multi-isoform locus panels", fontsize=13, weight="bold")
    save_figure(fig, OUT, "p02_isoforms")


def plot_regions() -> None:
    fig, axes = make_grid(1, figsize=(8, 2.4), constrained_layout=True)
    axis = axes[0, 0]
    intervals = [
        GenomicInterval("chr1", 120, 260), GenomicInterval("chr1", 290, 410),
        GenomicInterval("chr1", 670, 820),
    ]
    for row, index_row in enumerate(pack_intervals([(item.start, item.end) for item in intervals])):
        for index in index_row:
            item = intervals[index]
            axis.add_patch(
                plt.Rectangle((item.start, row), item.length, 0.7, color="#59a14f", alpha=0.8)
            )
    axis.set_xlim(0, 1000)
    axis.set_ylim(-0.2, 1.9)
    axis.set_yticks([0.35, 1.35], ["enhancers", "peaks"])
    axis.axhline(1.35, color="#59636e", linewidth=0.8)
    style_axis(axis, "Aligned genomic interval evidence")
    save_figure(fig, OUT, "p03_regions")


def plot_ratio() -> None:
    ratio = np.log2((SIGNAL_B + 0.05) / (SIGNAL_A + 0.05))
    fig, axes_grid = make_grid(2, figsize=(8, 4.2), sharex=True, constrained_layout=True)
    axes = axes_grid[:, 0]
    plot_signal(axes[0], X, SIGNAL_A, color="#377eb8", fill=False, linewidth=1.0, label="WT")
    plot_signal(axes[0], X, SIGNAL_B, color="#e41a1c", fill=False, linewidth=1.0, label="MUT")
    axes[0].legend(frameon=False, ncol=2, fontsize=8)
    axes[1].fill_between(X, 0, ratio, where=ratio >= 0, color="#d95f02", alpha=0.75)
    axes[1].fill_between(X, 0, ratio, where=ratio < 0, color="#1b9e77", alpha=0.75)
    axes[1].axhline(0, color="#404040", linewidth=0.8)
    style_axis(axes[0], "Matched signal comparison")
    style_axis(axes[1], "Explicit log2 ratio · pseudocount 0.05")
    axes[-1].set_xlabel("synthetic chromosome 1 coordinate")
    save_figure(fig, OUT, "p04_ratio")


def plot_highlights() -> None:
    fig, axes = make_grid(1, figsize=(8, 2.8), constrained_layout=True)
    axis = axes[0, 0]
    plot_signal(axis, X, SIGNAL_A, color="#2f5597", fill=False, linewidth=1.5)
    for start, end, label in [(240, 340, "candidate enhancer"), (660, 750, "mutant-enriched peak")]:
        axis.axvspan(start, end, color="#f2c14e", alpha=0.3)
        axis.text((start + end) / 2, 1.18, label, ha="center", fontsize=8)
    style_axis(axis, "Signal with non-destructive evidence highlights")
    axis.set_xlabel("synthetic chromosome 1 coordinate")
    save_figure(fig, OUT, "p05_highlights")


def plot_zoom() -> None:
    """Two-column parent/child zoom panels with complete genomic context."""
    fig, axes = make_grid(
        12, 2, figsize=(13, 12),
        height_ratios=[1, 1, 1, 1, 0.45, 0.35, 1, 1, 1, 1, 0.45, 0.35],
        constrained_layout=True,
    )
    for column, focus in enumerate([(250, 390), (610, 760)]):
        parent_x = np.linspace(0, 1000, 220)
        for track in range(4):
            axis = axes[track, column]
            values = (track + 1) * (
                np.exp(-((parent_x - (focus[0] + 30)) / 55) ** 2)
                + 0.4 * np.exp(-((parent_x - 700) / 45) ** 2)
            )
            plot_signal(axis, parent_x, values,
                        color="#c7c7c7" if track < 2 else "#555555",
                        alpha=1.0, linewidth=0)
            axis.axvspan(*focus, color="#f2c14e", alpha=0.22)
            axis.set_xlim(0, 1000)
            axis.set_yticks([])
            axis.spines[["top", "right", "left"]].set_visible(False)
            if track == 0:
                axis.set_title(
                    f"Parent locus {column + 1} | chr1:0-1,000",
                    loc="left", fontsize=10, weight="bold",
                )
            if track == 3:
                axis.set_xlabel("parent coordinate", fontsize=8)
        draw_gene_model(axes[4, column], reverse=column % 2 == 1, labels=False)
        connector = axes[5, column]
        connector.fill(
            [0.0, 1.0, 0.78, 0.22], [0.9, 0.9, 0.05, 0.05],
            color="#8c8c8c", alpha=0.65,
        )
        connector.plot([0.0, 0.22], [0.9, 0.05], color="#666666", linewidth=0.8)
        connector.plot([1.0, 0.78], [0.9, 0.05], color="#666666", linewidth=0.8)
        connector.set_xlim(0, 1)
        connector.set_ylim(0, 1)
        connector.set_axis_off()
        zoom_x = np.linspace(*focus, 180)
        for track in range(4):
            axis = axes[6 + track, column]
            values = (track + 1) * np.exp(-((zoom_x - np.mean(focus)) / 18) ** 2)
            plot_signal(axis, zoom_x, values,
                        color="#c7c7c7" if track < 2 else "#555555",
                        alpha=1.0, linewidth=0)
            axis.set_xlim(*focus)
            axis.set_yticks([])
            axis.spines[["top", "right", "left"]].set_visible(False)
            if track == 0:
                axis.set_title("Child zoom", loc="left", fontsize=9)
            if track == 3:
                axis.set_xlabel("zoom coordinate", fontsize=8)
        draw_gene_model(axes[10, column], reverse=column % 2 == 1, labels=False)
        draw_cytoband(axes[11, column], marker=0.5)
    fig.suptitle("Independent parent-to-child zoom composition", fontsize=13, weight="bold")
    save_figure(fig, OUT, "p06_zoom")


if __name__ == "__main__":
    plot_faceted()
    plot_isoforms()
    plot_regions()
    plot_ratio()
    plot_highlights()
    plot_zoom()
