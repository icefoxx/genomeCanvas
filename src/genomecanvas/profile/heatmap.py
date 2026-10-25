"""Sorted region heatmaps from profile matrices."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np

from ..render.plotting import save_figure
from .extract import ProfileMatrix
from .summarize import profile_summarize

SortBy = Literal["mean", "median"]


def profile_heatmap(
    matrices: list[ProfileMatrix],
    *,
    output: str | Path | None = None,
    sort_by: SortBy = "mean",
    top_profile: bool = True,
    cmap: str = "Blues",
    zmax: float | None = None,
) -> Path | None:
    """Draw one heatmap column per sample, optionally with a mean profile on top."""
    if not matrices:
        raise ValueError("profile_heatmap requires at least one matrix")
    order = _sort_index(matrices[0].values, sort_by)
    n_cols = len(matrices)
    height_ratios = [1.2, 4.0] if top_profile else [4.0]
    rows = 2 if top_profile else 1
    figure, axes = plt.subplots(
        rows,
        n_cols,
        figsize=(3.2 * n_cols, 6.0 if top_profile else 4.8),
        squeeze=False,
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.08, "wspace": 0.12},
    )
    summaries = profile_summarize(matrices) if top_profile else {}
    up, down = matrices[0].up, matrices[0].down
    for column, matrix in enumerate(matrices):
        ordered = matrix.values[order] if matrix.values.size else matrix.values
        vmax = zmax if zmax is not None else float(np.nanpercentile(ordered, 99) or 1.0)
        image_axis = axes[-1, column]
        image_axis.imshow(
            ordered,
            aspect="auto",
            cmap=cmap,
            vmin=0.0,
            vmax=max(vmax, 1e-9),
            interpolation="bilinear",
            extent=(-up, down, ordered.shape[0], 0),
        )
        image_axis.axvline(0, color="#888888", linewidth=0.7, linestyle="--")
        image_axis.set_xlabel("distance (bp)")
        image_axis.set_title(matrix.sample, fontsize=10)
        image_axis.set_yticks([])
        if column == 0:
            image_axis.set_ylabel("regions")
        else:
            image_axis.tick_params(labelleft=False)
        if top_profile:
            axis = axes[0, column]
            curve = summaries[matrix.sample]
            x = np.linspace(-up, down, curve.size) if curve.size else np.array([0.0])
            axis.plot(x, curve, color="#1f4e79", linewidth=1.4)
            axis.set_xlim(-up, down)
            axis.spines[["top", "right", "bottom"]].set_visible(False)
            axis.tick_params(labelbottom=False, bottom=False)
            if column:
                axis.tick_params(labelleft=False)
    if output is None:
        return None
    path = Path(output)
    save_figure(figure, path.parent, path.stem, close=True)
    return path.with_suffix(".png")


def _sort_index(values: np.ndarray, sort_by: SortBy) -> np.ndarray:  # type: ignore[type-arg]
    if values.size == 0:
        return np.array([], dtype=int)
    scores = np.nanmean(values, axis=1) if sort_by == "mean" else np.nanmedian(values, axis=1)
    return np.argsort(-scores)
