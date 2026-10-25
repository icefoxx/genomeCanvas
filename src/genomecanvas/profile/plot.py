"""Draw summarized profile line plots."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..render.plotting import save_figure


def profile_plot(
    summary: dict[str, np.ndarray],  # type: ignore[type-arg]
    *,
    up: int,
    down: int,
    output: str | Path | None = None,
    colors: dict[str, str] | None = None,
    xlabel: str = "distance (bp)",
    ylabel: str = "signal",
    title: str | None = None,
    labels: dict[str, str] | None = None,
    legend: bool = True,
    ln_v: bool = True,
) -> Path | None:
    """Plot one summarized profile curve per sample or condition."""
    figure, axis = plt.subplots(figsize=(6.5, 3.6), facecolor="white")
    palette = colors or {}
    default = ["#2f4f4f", "#d35400", "#2980b9", "#27ae60", "#8e44ad", "#c0392b"]
    rename = labels or {}
    for index, (name, values) in enumerate(summary.items()):
        if values.size == 0:
            continue
        x = np.linspace(-up, down, values.size)
        color = palette.get(name, default[index % len(default)])
        axis.plot(x, values, color=color, linewidth=1.6, label=rename.get(name, name))
    if ln_v:
        axis.axvline(0, color="#888888", linewidth=0.8, linestyle="--")
    axis.set_xlabel(xlabel, fontsize=11)
    axis.set_ylabel(ylabel, fontsize=11)
    if title:
        axis.set_title(title, fontsize=12, pad=8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(labelsize=9)
    if legend:
        axis.legend(frameon=False, fontsize=9)
    figure.tight_layout()
    if output is None:
        return None
    path = Path(output)
    save_figure(figure, path.parent, path.stem, close=True)
    return path.with_suffix(".png")
