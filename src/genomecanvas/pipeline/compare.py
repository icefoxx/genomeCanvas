"""Peak-center Mann-Whitney tests between BED groups × samples."""

from __future__ import annotations

import math
from itertools import combinations
from pathlib import Path

import numpy as np

from ..render.plotting import save_figure
from .matrix import EnrichmentMatrix


def peak_center_comparisons(
    matrix: EnrichmentMatrix,
    output_prefix: str | Path,
    *,
    flank_bins: int = 5,
    significance: float = 0.05,
    cmap: str = "viridis_r",
) -> list[Path]:
    """Write per-meta-group p-value tables and -log10(p) heatmaps."""
    if len(matrix.group_labels) < 2:
        return []
    written: list[Path] = []
    center = matrix.bins_per_sample // 2
    lo = max(0, center - flank_bins)
    hi = min(matrix.bins_per_sample, center + flank_bins + 1)
    meta_to_samples: dict[str, list[int]] = {}
    for index, meta in enumerate(matrix.meta_groups):
        meta_to_samples.setdefault(meta, []).append(index)
    for meta, sample_indices in meta_to_samples.items():
        keys: list[str] = []
        values: dict[str, np.ndarray] = {}  # type: ignore[type-arg]
        for group_i, group in enumerate(matrix.group_labels):
            row_start = matrix.group_boundaries[group_i]
            row_end = matrix.group_boundaries[group_i + 1]
            for sample_i in sample_indices:
                key = f"BED_{group}._.sample_{matrix.sample_labels[sample_i]}"
                keys.append(key)
                block = matrix.sample_block(sample_i)[row_start:row_end, lo:hi]
                values[key] = np.nanmean(block, axis=1) if block.size else np.array([])
        if len(keys) < 2:
            continue
        n = len(keys)
        pvals = np.ones((n, n), dtype=float)
        for i, j in combinations(range(n), 2):
            left, right = values[keys[i]], values[keys[j]]
            pvals[i, j] = pvals[j, i] = _mannwhitney(left, right)
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in meta)
        prefix = Path(f"{output_prefix}.peakCenter.{safe}")
        tsv = prefix.with_suffix(".pvalues.tsv")
        tsv.write_text(
            "\t".join(["group", *keys]) + "\n"
            + "\n".join(
                keys[row] + "\t" + "\t".join(f"{pvals[row, column]:.6g}" for column in range(n))
                for row in range(n)
            )
            + "\n",
            encoding="utf-8",
        )
        neglog = -np.log10(np.clip(pvals, 1e-300, 1.0))
        np.fill_diagonal(neglog, 0.0)
        import matplotlib.pyplot as plt

        figure, axis = plt.subplots(figsize=(max(6.0, 0.45 * n + 3), max(5.5, 0.45 * n + 2.5)))
        image = axis.imshow(neglog, cmap=cmap)
        axis.set_xticks(range(n), keys, rotation=45, ha="right", fontsize=6)
        axis.set_yticks(range(n), keys, fontsize=6)
        figure.colorbar(
            image, ax=axis, fraction=0.046,
            label=f"-log10(p); threshold={significance}",
        )
        save_figure(figure, prefix.parent, prefix.name + ".clustermap", close=True)
        written.extend([tsv, prefix.parent / f"{prefix.name}.clustermap.png"])
    return written


def _mannwhitney(left: np.ndarray, right: np.ndarray) -> float:  # type: ignore[type-arg]
    left = left[np.isfinite(left)]
    right = right[np.isfinite(right)]
    if left.size == 0 or right.size == 0:
        return 1.0
    try:
        from scipy.stats import mannwhitneyu

        return float(mannwhitneyu(left, right, alternative="two-sided").pvalue)
    except ImportError:
        return _mw_numpy(left, right)


def _mw_numpy(left: np.ndarray, right: np.ndarray) -> float:  # type: ignore[type-arg]
    combined = np.concatenate([left, right])
    ranks = combined.argsort().argsort().astype(float) + 1.0
    r1 = ranks[: left.size].sum()
    n1, n2 = left.size, right.size
    u1 = r1 - n1 * (n1 + 1) / 2
    mean = n1 * n2 / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    if sd == 0:
        return 1.0
    z = abs(u1 - mean) / sd
    return float(max(0.0, min(1.0, 2 * (1 - 0.5 * (1 + math.erf(z / math.sqrt(2)))))))
