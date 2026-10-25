"""Volcano plot for differential peak results."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..render.plotting import save_figure
from .diffpeak import DiffPeakResult


def volcano_plot(
    result: DiffPeakResult,
    *,
    output: str | Path | None = None,
    fdr: float = 0.1,
    up_color: str = "#d35400",
    down_color: str = "#1abc9c",
) -> Path | None:
    """Scatter log2 fold-change against -log10 p-value."""
    figure, axis = plt.subplots(figsize=(5.4, 4.4))
    y = -np.log10(np.clip(result.p_value, 1e-300, 1.0))
    significant = result.fdr < fdr
    axis.scatter(
        result.log2_fold_change[~significant],
        y[~significant],
        s=12,
        color="#bbbbbb",
        alpha=0.7,
        linewidths=0,
    )
    up = significant & (result.log2_fold_change > 0)
    down = significant & (result.log2_fold_change < 0)
    axis.scatter(result.log2_fold_change[up], y[up], s=16, color=up_color, linewidths=0)
    axis.scatter(result.log2_fold_change[down], y[down], s=16, color=down_color, linewidths=0)
    axis.set_xlabel("log2 fold change")
    axis.set_ylabel("-log10(p)")
    axis.set_title(f"{result.numerator} / {result.denominator}", fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    if output is None:
        return None
    path = Path(output)
    save_figure(figure, path.parent, path.stem, close=True)
    return path.with_suffix(".png")
