"""Principal component scores for region summaries."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..render.plotting import save_figure
from .summary import RegionSummary


@dataclass(frozen=True, slots=True)
class PCAResult:
    """Sample coordinates and variance explained."""

    samples: tuple[str, ...]
    coordinates: np.ndarray  # type: ignore[type-arg]
    variance_explained: np.ndarray  # type: ignore[type-arg]


def pca_scores(
    summary: RegionSummary, *, top: int = 1000, log2: bool = False, n_components: int = 2
) -> PCAResult:
    """Compute PCA on the most variable regions (samples as observations)."""
    values = np.asarray(summary.values, dtype=float)
    if log2:
        values = np.log2(values + 0.1)
    if values.shape[0] > top:
        order = np.argsort(np.nanstd(values, axis=1))[::-1][:top]
        values = values[order]
    table = np.nan_to_num(values.T, nan=0.0)
    centered = table - table.mean(axis=0, keepdims=True)
    _u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    coordinates = centered @ vt[: max(2, n_components)].T
    variance = singular ** 2
    explained = variance / variance.sum() if variance.sum() else variance
    return PCAResult(summary.samples, coordinates, explained)


def pca_plot(
    result: PCAResult,
    *,
    output: str | Path | None = None,
    colors: dict[str, str] | None = None,
    point_size: float = 48,
) -> Path | None:
    """Scatter the first two principal components."""
    figure, axis = plt.subplots(figsize=(5.2, 4.4))
    palette = colors or {}
    default = "#2c3e50"
    for index, sample in enumerate(result.samples):
        axis.scatter(
            result.coordinates[index, 0],
            result.coordinates[index, 1],
            color=palette.get(sample, default),
            s=point_size,
        )
        axis.annotate(
            sample, result.coordinates[index],
            textcoords="offset points", xytext=(5, 4), fontsize=8,
        )
    pc1 = result.variance_explained[0] if result.variance_explained.size else 0.0
    pc2 = result.variance_explained[1] if result.variance_explained.size > 1 else 0.0
    axis.set_xlabel(f"PC1 ({pc1:.2f})")
    axis.set_ylabel(f"PC2 ({pc2:.2f})")
    axis.axhline(0, color="#dddddd", linewidth=0.8)
    axis.axvline(0, color="#dddddd", linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    if output is None:
        return None
    path = Path(output)
    save_figure(figure, path.parent, path.stem, close=True)
    return path.with_suffix(".png")
