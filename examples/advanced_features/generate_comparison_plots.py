"""Real-data comparison-aware transformation figures."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))

from examples.reference_comparison.generate import DATA, SIGNALS
from examples.reference_comparison.generate_p02_p30 import LOCI7
from genomecanvas.io.bigwig import read_bigwig
from genomecanvas.model.background_spec import BackgroundSpec, TransformationSpec
from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.model.track import SignalSlice
from genomecanvas.provenance.manifest import build_manifest, write_manifest
from genomecanvas.render.plotting import make_grid, plot_signal, save_figure, style_axis
from genomecanvas.tracks.difference import compare_signal
from genomecanvas.transform.background import estimate_background, transform_values

OUT = Path(__file__).parent
PSEUDOCOUNT = 0.01


def slice_for(path: Path, chrom: str, start: int, end: int) -> SignalSlice:
    """Read a BigWig interval into a normalized SignalSlice."""
    values = np.asarray(read_bigwig(str(path), chrom, start, end), dtype=float)
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    starts = np.linspace(start, end, len(values), endpoint=False, dtype=int)
    ends = np.linspace(start, end, len(values) + 1, dtype=int)[1:]
    return SignalSlice(chrom, start, end, starts, ends, values, "signal_units")


def comparison_plot() -> None:
    """Plot parent tracks and an explicit B-A visual difference."""
    figure, axes = make_grid(3, 4, figsize=(13, 7))
    for column, (chrom, start, end) in enumerate(LOCI7):
        parent_a = slice_for(DATA / SIGNALS[0][1], chrom, start, end)
        parent_b = slice_for(DATA / SIGNALS[1][1], chrom, start, end)
        difference = compare_signal(parent_a, parent_b, operation="difference")
        for row, signal, label, color in [
            (0, parent_a, "Input rep1 · condition A", "#377eb8"),
            (1, parent_b, "Input rep2 · condition B", "#e41a1c"),
            (2, difference, "visual difference · B − A", "#984ea3"),
        ]:
            axis = axes[row, column]
            x = np.linspace(start, end, len(signal.values), endpoint=False)
            plot_signal(axis, x, signal.values, color=color, alpha=0.55)
            axis.set_xlim(start, end)
            axis.set_yticks([])
            style_axis(axis, hide_y=True)
            if column == 0:
                axis.set_ylabel(label, fontsize=8)
        axes[0, column].set_title(f"{chrom}:{start:,}-{end:,}", fontsize=8, weight="bold")
    figure.suptitle("Comparison-aware derived track: explicit B − A", weight="bold")
    figure.subplots_adjust(hspace=0.04, wspace=0.08, top=0.9)
    save_figure(figure, OUT, "p31_comparison_difference")


def ratio_plot() -> None:
    """Plot explicit ratio and log2 ratio with a declared pseudocount."""
    figure, axes = make_grid(2, 4, figsize=(13, 5))
    for column, (chrom, start, end) in enumerate(LOCI7):
        a = slice_for(DATA / SIGNALS[0][1], chrom, start, end)
        b = slice_for(DATA / SIGNALS[1][1], chrom, start, end)
        ratio = compare_signal(a, b, operation="ratio", pseudocount=PSEUDOCOUNT)
        log2 = compare_signal(a, b, operation="log2_ratio", pseudocount=PSEUDOCOUNT)
        for row, signal, label in [(0, ratio, "B/A"), (1, log2, "log2(B/A)")]:
            axis = axes[row, column]
            x = np.linspace(start, end, len(signal.values), endpoint=False)
            axis.plot(x, signal.values, color="#e66101", linewidth=0.7)
            axis.axhline(1 if row == 0 else 0, color="#555555", linewidth=0.6)
            axis.set_xlim(start, end)
            axis.set_yticks([])
            style_axis(axis, hide_y=True)
            if column == 0:
                axis.set_ylabel(label, fontsize=8)
        axes[0, column].set_title(f"{chrom}:{start:,}-{end:,}", fontsize=8)
    figure.suptitle(f"Explicit ratio transformations · pseudocount {PSEUDOCOUNT}", weight="bold")
    figure.subplots_adjust(hspace=0.05, wspace=0.08, top=0.88)
    save_figure(figure, OUT, "p32_comparison_ratios")


def local_background_plot() -> None:
    """Plot raw signal against local-background-adjusted enrichment."""
    figure, axes = make_grid(2, 4, figsize=(13, 5))
    background_spec = BackgroundSpec(
        upstream=10_000, downstream=10_000, estimator="median"
    )
    for column, (chrom, start, end) in enumerate(LOCI7):
        signal = slice_for(DATA / SIGNALS[2][1], chrom, start, end)
        left = slice_for(DATA / SIGNALS[2][1], chrom, max(0, start - 10_000), start)
        right = slice_for(DATA / SIGNALS[2][1], chrom, end, end + 10_000)
        estimate = estimate_background(left.values, right.values, background_spec)
        adjusted = transform_values(signal.values, estimate.value, "log2_ratio", PSEUDOCOUNT)
        x = np.linspace(start, end, len(signal.values), endpoint=False)
        plot_signal(axes[0, column], x, signal.values, color="#555555", alpha=0.55)
        plot_signal(axes[1, column], x, adjusted, color="#7b3294", alpha=0.65)
        axes[1, column].axhline(0, color="#555555", linewidth=0.6)
        for row in range(2):
            axes[row, column].set_xlim(start, end)
            axes[row, column].set_yticks([])
            style_axis(axes[row, column], hide_y=True)
        axes[0, column].set_title(f"background={estimate.value:.2f}", fontsize=8)
    axes[0, 0].set_ylabel("raw IP rep1", fontsize=8)
    axes[1, 0].set_ylabel("log2(local enrichment)", fontsize=8)
    figure.suptitle("Local-context transformation · not sample normalization", weight="bold")
    figure.subplots_adjust(hspace=0.04, wspace=0.08, top=0.88)
    save_figure(figure, OUT, "p33_local_background_comparison")


def provenance_example() -> None:
    """Write a FigureSpec and manifest describing the comparison plot."""
    spec = FigureSpec.model_validate({
        "genome": {"assembly": "GRCm38"},
        "region": {"chrom": "chr7", "start": LOCI7[0][1], "end": LOCI7[0][2]},
        "tracks": [],
        "transformations": [
            TransformationSpec(
                id="local_enrichment", operation="log2_ratio", input="IP_rep1",
                background=BackgroundSpec(upstream=10_000, downstream=10_000),
                pseudocount=PSEUDOCOUNT,
            )
        ],
    })
    (OUT / "p33_local_background.spec.yaml").write_text(spec.to_yaml(), encoding="utf-8")
    manifest = build_manifest(spec)
    write_manifest(manifest, OUT / "p33_local_background.manifest.json")
    (OUT / "p33_local_background.summary.json").write_text(
        json.dumps({"figure_spec": spec.to_dict(), "manifest": manifest}, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    comparison_plot()
    ratio_plot()
    local_background_plot()
    provenance_example()
