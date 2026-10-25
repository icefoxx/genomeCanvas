"""Two-group differential peak testing on region summaries."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .summary import RegionSummary


@dataclass(frozen=True, slots=True)
class DiffPeakResult:
    """Per-region Welch t-test with Benjamini-Hochberg FDR."""

    regions: tuple[str, ...]
    log2_fold_change: np.ndarray  # type: ignore[type-arg]
    p_value: np.ndarray  # type: ignore[type-arg]
    fdr: np.ndarray  # type: ignore[type-arg]
    numerator: str
    denominator: str


def diffpeak(
    summary: RegionSummary,
    conditions: dict[str, str],
    *,
    numerator: str,
    denominator: str,
    log2: bool = True,
) -> DiffPeakResult:
    """Compare two conditions with a Welch t-test on log2 intensities."""
    values = np.asarray(summary.values, dtype=float)
    if log2:
        values = np.log2(values + 0.1)
    num_idx = [
        index for index, sample in enumerate(summary.samples)
        if conditions[sample] == numerator
    ]
    den_idx = [
        index for index, sample in enumerate(summary.samples)
        if conditions[sample] == denominator
    ]
    if len(num_idx) < 1 or len(den_idx) < 1:
        raise ValueError("each condition must contain at least one sample")
    left = values[:, den_idx]
    right = values[:, num_idx]
    mean_left = np.nanmean(left, axis=1)
    mean_right = np.nanmean(right, axis=1)
    fold = mean_right - mean_left
    p_values = np.array(
        [_welch_p(right[row], left[row]) for row in range(values.shape[0])],
        dtype=float,
    )
    return DiffPeakResult(
        regions=summary.regions,
        log2_fold_change=fold,
        p_value=p_values,
        fdr=_benjamini_hochberg(p_values),
        numerator=numerator,
        denominator=denominator,
    )


def _welch_p(a: np.ndarray, b: np.ndarray) -> float:  # type: ignore[type-arg]
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return 1.0
    mean_a, mean_b = float(a.mean()), float(b.mean())
    var_a = float(a.var(ddof=1)) if a.size > 1 else 0.0
    var_b = float(b.var(ddof=1)) if b.size > 1 else 0.0
    denom = np.sqrt(var_a / max(a.size, 1) + var_b / max(b.size, 1))
    if denom == 0:
        return 1.0 if mean_a == mean_b else 0.0
    t_stat = abs(mean_a - mean_b) / denom
    p_value = 2.0 * (1.0 - _phi(t_stat))
    return float(min(1.0, max(0.0, p_value)))


def _phi(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:  # type: ignore[type-arg]
    n = p_values.size
    if n == 0:
        return p_values
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(n, dtype=float)
    result[order] = np.clip(adjusted, 0, 1)
    return result
