"""Per-meta-group display renormalization of enrichment matrices."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from .matrix import EnrichmentMatrix


def renormalize(
    matrix: EnrichmentMatrix, *, approach: str = "profileMinMax"
) -> EnrichmentMatrix:
    """Return a copy scaled by ``profileMinMax`` or ``qnorm`` within meta groups."""
    values = np.array(matrix.values, dtype=np.float64, copy=True)
    groups: dict[str, list[int]] = defaultdict(list)
    for index, name in enumerate(matrix.meta_groups):
        groups[name].append(index)
    if approach == "qnorm":
        values = _quantile_normalize(values, matrix, groups)
        approach = "profileMinMax"
    if approach != "profileMinMax":
        raise ValueError("renorm approach must be profileMinMax or qnorm")
    profiles = {
        label: np.nanmean(matrix.sample_block(index), axis=0)
        for index, label in enumerate(matrix.sample_labels)
    }
    for _group, indices in groups.items():
        if not indices:
            continue
        stacked = np.vstack([profiles[matrix.sample_labels[index]] for index in indices])
        low, high = float(np.nanmin(stacked)), float(np.nanmax(stacked))
        scale = high - low if high > low else 1.0
        for index in indices:
            start, end = matrix.sample_boundaries[index], matrix.sample_boundaries[index + 1]
            values[:, start:end] = (values[:, start:end] - low) / scale
    return EnrichmentMatrix(
        values=np.asarray(values, dtype=np.float32),
        coords=matrix.coords,
        sample_labels=matrix.sample_labels,
        sample_paths=matrix.sample_paths,
        sample_boundaries=matrix.sample_boundaries,
        group_labels=matrix.group_labels,
        group_boundaries=matrix.group_boundaries,
        upstream=matrix.upstream,
        downstream=matrix.downstream,
        body=matrix.body,
        bin_size=matrix.bin_size,
        ref_point=matrix.ref_point,
        bins_per_sample=matrix.bins_per_sample,
        meta_groups=matrix.meta_groups,
        conditions=matrix.conditions,
        colormaps=matrix.colormaps,
        profile_colors=matrix.profile_colors,
    )


def _quantile_normalize(
    values: np.ndarray,  # type: ignore[type-arg]
    matrix: EnrichmentMatrix,
    groups: dict[str, list[int]],
) -> np.ndarray:  # type: ignore[type-arg]
    out = np.array(values, dtype=np.float64, copy=True)
    n_regions = values.shape[0]
    for indices in groups.values():
        if len(indices) <= 1:
            continue
        flats = []
        spans = []
        for index in indices:
            start, end = matrix.sample_boundaries[index], matrix.sample_boundaries[index + 1]
            block = values[:, start:end]
            flats.append(block.reshape(-1))
            spans.append((start, end))
        stacked = np.column_stack(flats)
        sorted_cols = np.sort(stacked, axis=0)
        rank_means = np.mean(sorted_cols, axis=1)
        ranks = np.argsort(np.argsort(stacked, axis=0), axis=0)
        for column, (start, end) in enumerate(spans):
            restored = rank_means[ranks[:, column]].reshape(n_regions, end - start)
            out[:, start:end] = restored
    return out
