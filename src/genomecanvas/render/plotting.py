"""Shared Matplotlib helpers for GenomeCanvas example figures."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt  # type: ignore[import-untyped]
import numpy as np


def make_grid(
    rows: int,
    columns: int = 1,
    *,
    figsize: tuple[float, float],
    sharex: bool | str = False,
    height_ratios: Sequence[float] | None = None,
    constrained_layout: bool = False,
    gridspec_kw: dict[str, Any] | None = None,
) -> tuple[plt.Figure, np.ndarray[object, object]]:
    """Create every example figure with a predictable two-dimensional axes array."""
    grid_options = dict(gridspec_kw or {})
    if height_ratios is not None:
        grid_options["height_ratios"] = height_ratios
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=figsize,
        sharex=sharex,
        squeeze=False,
        constrained_layout=constrained_layout,
        gridspec_kw=grid_options or None,
    )
    return figure, axes


def save_figure(
    figure: plt.Figure,
    output_dir: Path,
    name: str,
    *,
    dpi: int = 180,
    close: bool = True,
) -> None:
    """Write the standard PNG/PDF pair and release the Matplotlib figure."""
    figure.savefig(output_dir / f"{name}.png", dpi=dpi, bbox_inches="tight")
    figure.savefig(output_dir / f"{name}.pdf", bbox_inches="tight")
    if close:
        plt.close(figure)


def style_axis(
    axis: plt.Axes,
    title: str | None = None,
    *,
    grid: bool = False,
    hide_y: bool = False,
    hide_x: bool = False,
) -> None:
    """Apply the common clean frame used by example signal panels."""
    if title:
        axis.set_title(title, loc="left", fontsize=10, weight="bold")
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0, labelsize=8)
    if grid:
        axis.grid(axis="y", color="#e8edf2", linewidth=0.6)
    if hide_y:
        axis.set_yticks([])
        axis.tick_params(axis="y", left=False, labelleft=False)
    if hide_x:
        axis.tick_params(axis="x", bottom=False, labelbottom=False)


def plot_signal(
    axis: plt.Axes,
    x: np.ndarray,
    values: np.ndarray,
    *,
    color: str,
    fill: bool = True,
    alpha: float = 0.55,
    linewidth: float = 1.0,
    baseline: float | None = None,
    label: str | None = None,
) -> None:
    """Draw a signal consistently as a line, optionally with an area fill."""
    if fill:
        bottom = 0.0 if baseline is None else baseline
        axis.fill_between(x, bottom, values, color=color, alpha=alpha, linewidth=0)
    axis.plot(x, values, color=color, linewidth=linewidth, label=label)
