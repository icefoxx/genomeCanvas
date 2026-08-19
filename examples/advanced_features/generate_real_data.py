"""Real-data plots for local-background-adjusted visualization."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))

from examples.reference_comparison.generate import DATA, SIGNALS
from examples.reference_comparison.generate_p02_p30 import LOCI7, LOCI14
from genomecanvas.io.bigwig import read_bigwig
from genomecanvas.model.background_spec import BackgroundSpec
from genomecanvas.render.plotting import make_grid, plot_signal, save_figure, style_axis
from genomecanvas.transform.background import estimate_background, transform_values

OUT = Path(__file__).parent
FLANK = 10_000
PSEUDOCOUNT = 0.01
TRACK = SIGNALS[2]


def values(path: Path, chrom: str, start: int, end: int) -> np.ndarray:
    """Read one requested BigWig interval through the GenomeCanvas adapter."""
    return np.nan_to_num(
        np.asarray(read_bigwig(str(path), chrom, start, end), dtype=float),
        nan=0.0, posinf=0.0, neginf=0.0,
    )


def adjusted(path: Path, chrom: str, start: int, end: int) -> tuple[np.ndarray, float]:
    """Return raw values and explicit median-flank log2 enrichment."""
    raw = values(path, chrom, start, end)
    left = values(path, chrom, max(0, start - FLANK), start)
    right = values(path, chrom, end, end + FLANK)
    estimate = estimate_background(
        left, right, BackgroundSpec(upstream=FLANK, downstream=FLANK, estimator="median")
    )
    transformed = transform_values(raw, estimate.value, "log2_ratio", PSEUDOCOUNT)
    return transformed, estimate.value


def render(name: str, loci: list[tuple[str, int, int]]) -> None:
    """Render raw and local-enrichment rows for several real loci."""
    path = DATA / TRACK[1]
    figure, axes = make_grid(2, len(loci), figsize=(13, 5))
    for column, (chrom, start, end) in enumerate(loci):
        raw = values(path, chrom, start, end)
        local, background = adjusted(path, chrom, start, end)
        x = np.linspace(start, end, len(raw), endpoint=False)
        plot_signal(axes[0, column], x, raw, color="#555555", alpha=1.0, linewidth=0)
        plot_signal(axes[1, column], x, local, color="#7b3294", alpha=1.0, linewidth=0)
        axes[0, column].set_title(f"{chrom}:{start:,}-{end:,}", fontsize=8, weight="bold")
        axes[0, column].text(
            0.98, 0.92, f"background={background:.2f}",
            transform=axes[0, column].transAxes, ha="right", fontsize=7,
        )
        for row in range(2):
            style_axis(axes[row, column], hide_x=True)
        axes[1, column].set_xlabel("shared locus", fontsize=8)
    axes[0, 0].set_ylabel("Raw signal", fontsize=9)
    axes[1, 0].set_ylabel("Log2(local enrichment)", fontsize=9)
    figure.suptitle(
        f"{name} · {TRACK[0]} · median flanks ±{FLANK:,} bp · pseudocount {PSEUDOCOUNT}",
        fontsize=11, weight="bold",
    )
    figure.subplots_adjust(hspace=0.08, wspace=0.12, top=0.84, left=0.06, right=0.98)
    save_figure(figure, OUT, f"{name}_local_background")


if __name__ == "__main__":
    render("p01", LOCI7)
    render("p02", LOCI14)
