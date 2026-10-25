"""Sample-level correlation, PCA, and pairwise scatter extras."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.pca import PCAResult, pca_plot
from ..analysis.summary import RegionSummary
from ..render.plotting import save_figure


def region_means_table(values: np.ndarray, sample_boundaries: tuple[int, ...]) -> np.ndarray:  # type: ignore[type-arg]
    """Reduce each sample's bins to one mean per region."""
    columns = []
    for start, end in zip(sample_boundaries[:-1], sample_boundaries[1:], strict=True):
        columns.append(np.nanmean(values[:, start:end], axis=1))
    return np.column_stack(columns) if columns else np.zeros((values.shape[0], 0))


def correlation_heatmap(
    table: np.ndarray,  # type: ignore[type-arg]
    labels: tuple[str, ...],
    output: str | Path,
    *,
    method: str = "pearson",
    cmap: str = "bwr",
    annotate: bool = True,
) -> Path:
    """Ward-ordered sample correlation heatmap."""
    if table.shape[1] < 2:
        raise ValueError("correlation heatmap needs at least two samples")
    if method == "spearman":
        ranks = np.vstack([_rankdata(table[:, column]) for column in range(table.shape[1])]).T
        corr = np.corrcoef(ranks, rowvar=False)
    else:
        corr = np.corrcoef(table, rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0)
    order = _ward_order(corr)
    ordered = corr[np.ix_(order, order)]
    names = [labels[index].replace("\n", " ") for index in order]
    width = max(5.0, 0.7 * len(names) + 2)
    height = max(4.5, 0.7 * len(names) + 1.5)
    figure, axis = plt.subplots(figsize=(width, height))
    image = axis.imshow(ordered, cmap=cmap, vmin=-1, vmax=1)
    axis.set_xticks(range(len(names)), names, rotation=45, ha="right", fontsize=8)
    axis.set_yticks(range(len(names)), names, fontsize=8)
    for row in range(len(names) if annotate else 0):
        for column in range(len(names)):
            axis.text(
                column, row, f"{ordered[row, column]:.2f}",
                ha="center", va="center", fontsize=7,
                color="white" if abs(ordered[row, column]) > 0.55 else "#1f2933",
            )
    figure.colorbar(image, ax=axis, fraction=0.046, label=f"{method} correlation")
    path = Path(output)
    save_figure(figure, path.parent, path.stem, close=True)
    path.with_suffix(".tsv").write_text(
        "\t".join(["sample", *names]) + "\n"
        + "\n".join(
            names[row]
            + "\t"
            + "\t".join(f"{ordered[row, col]:.6g}" for col in range(len(names)))
            for row in range(len(names))
        )
        + "\n",
        encoding="utf-8",
    )
    return path.with_suffix(".png")


def pair_grid(
    table: np.ndarray,  # type: ignore[type-arg]
    labels: tuple[str, ...],
    output: str | Path,
    *,
    max_rows: int = 5000,
    color: str = "#343d46",
    hist_color: str = "#dfe2e6",
) -> Path:
    """Lower-triangle scatter, diagonal histograms, upper-triangle correlations."""
    n = table.shape[1]
    if table.shape[0] > max_rows:
        rng = np.random.default_rng(7201)
        table = table[rng.choice(table.shape[0], max_rows, replace=False)]
    logged = np.log10(np.clip(table, 0, None) + 1.0)
    figure, axes = plt.subplots(n, n, figsize=(2.2 * n, 2.2 * n), squeeze=False)
    names = [label.replace("\n", " ") for label in labels]
    for row in range(n):
        for column in range(n):
            axis = axes[row, column]
            if row == column:
                axis.hist(logged[:, row], bins=20, color=hist_color, edgecolor="#65737e")
            elif row > column:
                axis.scatter(
                    logged[:, column], logged[:, row],
                    s=4, color=color, alpha=0.4, linewidths=0,
                )
            else:
                if logged.shape[0] > 2:
                    corr = np.corrcoef(logged[:, column], logged[:, row])[0, 1]
                else:
                    corr = float("nan")
                axis.text(
                    0.5, 0.5, f"r={corr:.2f}",
                    ha="center", va="center", transform=axis.transAxes,
                )
                axis.set_xticks([])
                axis.set_yticks([])
            if row == n - 1:
                axis.set_xlabel(names[column], fontsize=7)
            if column == 0:
                axis.set_ylabel(names[row], fontsize=7)
            axis.tick_params(labelsize=6)
    figure.suptitle(f"PairGrid log10(signal+1) for {logged.shape[0]} regions", fontsize=10)
    figure.tight_layout()
    path = Path(output)
    save_figure(figure, path.parent, path.stem, close=True)
    return path.with_suffix(".png")


def sample_pca(table: np.ndarray, labels: tuple[str, ...], output: str | Path) -> Path:  # type: ignore[type-arg]
    """PCA of samples using region-mean values."""
    summary = RegionSummary(
        regions=tuple(str(index) for index in range(table.shape[0])),
        samples=labels,
        chroms=tuple("." for _ in range(table.shape[0])),
        starts=tuple(0 for _ in range(table.shape[0])),
        ends=tuple(1 for _ in range(table.shape[0])),
        values=np.asarray(table, dtype=float),
    )
    from ..analysis.pca import pca_scores

    result: PCAResult = pca_scores(summary, log2=True)
    path = pca_plot(result, output=output)
    return Path(path) if path else Path(output)


def _rankdata(values: np.ndarray) -> np.ndarray:  # type: ignore[type-arg]
    order = np.argsort(values)
    ranks = np.empty(values.size, dtype=float)
    ranks[order] = np.arange(values.size)
    return ranks


def _ward_order(corr: np.ndarray) -> np.ndarray:  # type: ignore[type-arg]
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform
    except ImportError:
        return np.arange(corr.shape[0])
    distance = np.clip(1 - corr, 0, None)
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=False)
    tree = linkage(condensed, method="ward")
    return np.asarray(leaves_list(tree), dtype=int)
