"""Generate visual demonstrations of TASK-101–108 features."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from genomecanvas.model.background_spec import BackgroundSpec
from genomecanvas.model.scale import ScaleSpec
from genomecanvas.model.track import SignalSlice
from genomecanvas.render.plotting import make_grid, plot_signal, save_figure, style_axis
from genomecanvas.transform.automatic_focus import largest_peak
from genomecanvas.transform.background import estimate_background, transform_values

OUT = Path(__file__).parent
X = np.arange(0, 1000)
BASELINE = np.full(X.size, 10.0)
SIGNAL = BASELINE + 100 * np.exp(-((X - 330) / 45) ** 2) + 35 * np.exp(-((X - 760) / 70) ** 2)
LEFT = BASELINE.copy()
RIGHT = BASELINE.copy()
spec = BackgroundSpec(upstream=100, downstream=100, estimator="median")
estimate = estimate_background(LEFT, RIGHT, spec)


def plot_background_transforms() -> None:
    """Show raw, subtraction, ratio, and log2-ratio tracks."""
    subtract = transform_values(SIGNAL, estimate.value, "subtract")
    ratio = transform_values(SIGNAL, estimate.value, "ratio", 0.01)
    log2_ratio = transform_values(SIGNAL, estimate.value, "log2_ratio", 0.01)
    fig, axes_grid = make_grid(4, figsize=(9, 7), sharex=True)
    axes = axes_grid[:, 0]
    series = [
        (SIGNAL, "Raw signal", "#555555"),
        (subtract, "Local background-adjusted signal · subtract", "#377eb8"),
        (ratio, "Local enrichment · ratio", "#e66101"),
        (log2_ratio, "Log2(local enrichment) · pseudocount 0.01", "#984ea3"),
    ]
    for axis, (values, label, color) in zip(axes, series, strict=True):
        plot_signal(axis, X, values, color=color, alpha=0.35)
        axis.set_ylabel(label, fontsize=8)
        style_axis(axis)
    axes[-1].set_xlabel("synthetic genomic coordinate")
    fig.suptitle(f"Local background estimate: median = {estimate.value:.1f}", weight="bold")
    fig.subplots_adjust(hspace=0.05, top=0.92)
    save_figure(fig, OUT, "local_background_transforms")


def plot_scaling_and_focus() -> None:
    """Show fixed/robust display scales and automatic largest-peak focus."""
    starts = np.arange(0, 1000)
    ends = starts + 1
    slice_ = SignalSlice("chr1", 0, 1000, starts, ends, SIGNAL, "visual units")
    focus = largest_peak(slice_, 180)
    fixed_limits = ScaleSpec("fixed", ymin=0, ymax=160).resolve(SIGNAL)
    robust_limits = ScaleSpec("robust", quantile=0.95).resolve(SIGNAL)
    fig, axes_grid = make_grid(2, figsize=(9, 4.5), sharex=True)
    axes = axes_grid[:, 0]
    plot_signal(axes[0], X, SIGNAL, color="#555555", fill=False)
    axes[0].set_ylim(*fixed_limits)
    axes[0].axvspan(focus.start, focus.end, color="#fdae61", alpha=0.35)
    axes[0].set_ylabel("Fixed scale 0–160")
    axes[0].set_title("Automatic largest-peak highlight and fixed display scale", loc="left")
    plot_signal(axes[1], X, SIGNAL, color="#1b7837", fill=False)
    axes[1].set_ylim(*robust_limits)
    axes[1].axvspan(focus.start, focus.end, color="#fdae61", alpha=0.35)
    axes[1].set_ylabel("Robust 95% scale")
    axes[1].set_xlabel("synthetic genomic coordinate")
    for axis in axes:
        style_axis(axis)
    fig.tight_layout()
    save_figure(fig, OUT, "scaling_and_automatic_focus")


if __name__ == "__main__":
    plot_background_transforms()
    plot_scaling_and_focus()
