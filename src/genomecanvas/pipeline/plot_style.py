"""Publication-oriented plot chrome, modelled on SeqPlots label/key controls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np

PlotScale = Literal["linear", "log2", "zscore"]
DistanceUnit = Literal["auto", "bp", "kb", "Mb"]
LegendLoc = Literal[
    "best", "upper right", "upper left", "lower left", "lower right",
    "right", "center left", "center right", "lower center", "upper center", "center",
]


@dataclass
class HeatmapPlotStyle:
    """SeqPlots-like knobs for enrichment heatmaps with optional mean profiles."""

    title: str | None = None
    labels: Sequence[str] | None = None
    xlab: str | None = None
    ylab: str = "regions"
    colorbar_label: str = "signal"
    show_profile: bool = True
    show_colorbar: bool = True
    independent_scale: bool = True
    autoscale: bool = True
    zmin: float | None = None
    zmax: float | None = None
    plot_scale: PlotScale = "linear"
    ln_v: bool = True
    ln_h: bool = True
    interpolation: str = "bilinear"
    font_family: str = "DejaVu Sans"
    title_size: float = 12.0
    axis_label_size: float = 10.0
    tick_size: float = 8.0
    sample_label_size: float = 10.0
    legend_size: float = 8.0
    group_label_size: float = 8.0
    height: float = 8.0
    width_per_sample: float = 2.7
    profile_ratio: float = 0.28
    colorbar_ratio: float = 0.06
    linewidth: float = 1.5
    error_estimates: bool = True
    ref_label: str | None = None
    x_unit: DistanceUnit = "auto"
    show_row_numbers: bool = False
    show_group_labels: bool = True
    profile_ymin: float | None = None
    profile_ymax: float | None = None
    cmap: str | None = None
    dpi: int = 180


@dataclass
class ProfilePlotStyle:
    """SeqPlots-like knobs for average/profile line plots."""

    title: str | None = None
    labels: Sequence[str] | None = None
    xlab: str | None = None
    ylab: str = "mean signal"
    legend: bool = True
    legend_pos: LegendLoc = "best"
    error_estimates: bool = True
    ln_v: bool = True
    ln_h: float | None = None
    plot_scale: PlotScale = "linear"
    xlim: tuple[float, float] | None = None
    ylim: tuple[float, float] | None = None
    colors: Sequence[str] | None = None
    font_family: str = "DejaVu Sans"
    title_size: float = 12.0
    axis_label_size: float = 11.0
    tick_size: float = 9.0
    legend_size: float = 9.0
    linewidth: float = 1.8
    fill_alpha: float = 0.18
    x_unit: DistanceUnit = "auto"
    ref_label: str | None = None
    figsize: tuple[float, float] = (7.2, 4.2)
    dpi: int = 180


def scale_values(values: np.ndarray, mode: PlotScale) -> np.ndarray:  # type: ignore[type-arg]
    """Apply linear, log2, or z-score display scaling."""
    data = np.asarray(values, dtype=float)
    if mode == "log2":
        scaled = np.log2(np.clip(data, 0, None) + 1.0)
        scaled[~np.isfinite(scaled)] = np.nan
        return scaled
    if mode == "zscore":
        finite = data[np.isfinite(data)]
        if finite.size < 2:
            return data
        center = float(np.nanmean(finite))
        spread = float(np.nanstd(finite))
        if spread <= 0:
            return data - center
        return (data - center) / spread
    return data


def format_distance(value: float, unit: DistanceUnit = "auto") -> str:
    """Format a signed distance for an axis tick."""
    resolved = unit
    magnitude = abs(value)
    if unit == "auto":
        if magnitude >= 1_000_000:
            resolved = "Mb"
        elif magnitude >= 1000:
            resolved = "kb"
        else:
            resolved = "bp"
    if resolved == "Mb":
        number = value / 1_000_000
        text = f"{number:.2f}".rstrip("0").rstrip(".")
        return f"{text} Mb"
    if resolved == "kb":
        number = value / 1000
        if abs(number - round(number)) < 1e-9:
            return f"{int(round(number))} kb"
        text = f"{number:.2f}".rstrip("0").rstrip(".")
        return f"{text} kb"
    if abs(value - round(value)) < 1e-9:
        return f"{int(round(value))}"
    return f"{value:g}"


def format_genomic_position(position: float, span: int, unit: DistanceUnit = "auto") -> str:
    """Format an absolute genomic coordinate."""
    resolved = unit
    if unit == "auto":
        if span >= 5_000_000:
            resolved = "Mb"
        elif span >= 20_000:
            resolved = "kb"
        else:
            resolved = "bp"
    if resolved == "Mb":
        text = f"{position / 1_000_000:.2f}".rstrip("0").rstrip(".")
        return f"{text} Mb"
    if resolved == "kb":
        text = f"{position / 1000:.1f}".rstrip("0").rstrip(".")
        return f"{text} kb"
    return f"{int(round(position)):,}"


def mean_and_sem(block: np.ndarray) -> tuple[np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    """Return column mean and standard error, ignoring NaN."""
    mean = np.nanmean(block, axis=0)
    if block.shape[0] < 2:
        return mean, np.zeros_like(mean)
    std = np.nanstd(block, axis=0, ddof=1)
    count = np.maximum(np.sum(np.isfinite(block), axis=0), 1)
    return mean, std / np.sqrt(count)
