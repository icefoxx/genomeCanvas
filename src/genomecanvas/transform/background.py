"""Local genomic background estimation for visualization-only transforms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import numpy as np

from ..core.intervals import GenomicInterval
from ..model.background_spec import BackgroundSpec


@dataclass(frozen=True, slots=True)
class BackgroundEstimate:
    """Estimated background and diagnostics for a target region."""

    value: float
    left_value: float
    right_value: float
    estimator: str
    excluded_count: int = 0
    warnings: tuple[str, ...] = ()


def flanking_intervals(
    target: GenomicInterval, upstream: int, downstream: int
) -> tuple[GenomicInterval | None, GenomicInterval | None]:
    """Return clipped left/right flanks around a target interval."""
    if upstream < 0 or downstream < 0:
        raise ValueError("flank lengths must be non-negative")
    left = (
        GenomicInterval(target.chrom, max(0, target.start - upstream), target.start)
        if upstream else None
    )
    right = (
        GenomicInterval(target.chrom, target.end, target.end + downstream)
        if downstream else None
    )
    return left, right


def _finite(values: np.ndarray[Any, Any]) -> np.ndarray[Any, Any]:
    """Return finite values or raise a scientifically meaningful error."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError("local background is unavailable because it contains no finite values")
    return cast(np.ndarray[Any, Any], finite)


def estimate_values(
    values: np.ndarray[Any, Any], estimator: str = "median", trim_fraction: float = 0.1
) -> float:
    """Estimate a scalar background from finite signal values."""
    finite = _finite(values)
    if estimator == "mean":
        return float(np.mean(finite))
    if estimator == "median":
        return float(np.median(finite))
    if estimator == "trimmed_mean":
        if not 0 <= trim_fraction < 0.5:
            raise ValueError("trim_fraction must be in [0, 0.5)")
        count = int(np.floor(finite.size * trim_fraction))
        ordered = np.sort(finite)
        trimmed = ordered[count: finite.size - count] if count else ordered
        return float(np.mean(trimmed))
    raise ValueError("estimator must be mean, median, or trimmed_mean")


def estimate_background(
    left_values: np.ndarray[Any, Any],
    right_values: np.ndarray[Any, Any],
    spec: BackgroundSpec,
) -> BackgroundEstimate:
    """Estimate flanking background and emit asymmetry diagnostics."""
    if spec.strategy != "flanking":
        raise ValueError("estimate_background requires the flanking strategy")
    if spec.upstream is None or spec.downstream is None:
        raise ValueError("flanking background requires upstream and downstream lengths")
    left = estimate_values(left_values, spec.estimator, spec.trim_fraction)
    right = estimate_values(right_values, spec.estimator, spec.trim_fraction)
    value = estimate_values(
        np.concatenate([_finite(left_values), _finite(right_values)]),
        spec.estimator,
        spec.trim_fraction,
    )
    warnings: list[str] = []
    if spec.warning_threshold is not None and abs(left - right) > spec.warning_threshold:
        warnings.append("left and right flanking background estimates differ beyond threshold")
    if abs(value) <= np.finfo(float).eps:
        warnings.append("background is close to zero; ratio stability may be poor")
    if len(spec.exclusions) > 10:
        warnings.append("background contains many excluded intervals")
    return BackgroundEstimate(
        value, left, right, spec.estimator, len(spec.exclusions), tuple(warnings)
    )


def rolling_background(
    values: np.ndarray[Any, Any], window: int, estimator: str = "median"
) -> np.ndarray[Any, Any]:
    """Estimate a rolling local background without mutating source values."""
    if window < 1:
        raise ValueError("rolling background window must be positive")
    output = np.full(values.shape, np.nan, dtype=float)
    half = window // 2
    for index in range(len(values)):
        left, right = max(0, index - half), min(len(values), index + half + 1)
        neighborhood = np.concatenate([values[left:index], values[index + 1:right]])
        if neighborhood.size:
            output[index] = estimate_values(neighborhood, estimator)
    return output


def local_zscore(
    values: np.ndarray[Any, Any], background: float, scale: float
) -> np.ndarray[Any, Any]:
    """Compute a display-only local z-score with explicit scale."""
    if scale <= 0 or not np.isfinite(scale):
        raise ValueError("local z-score scale must be finite and positive")
    return cast(np.ndarray[Any, Any], (values - background) / scale)


def transform_values(
    signal_values: np.ndarray[Any, Any],
    background: float,
    operation: str,
    pseudocount: float | None = None,
) -> np.ndarray[Any, Any]:
    """Apply an explicit visualization transformation without mutating input."""
    if operation == "subtract":
        return signal_values - background
    if pseudocount is None:
        raise ValueError("ratio and log2_ratio require an explicit pseudocount")
    if pseudocount < 0:
        raise ValueError("pseudocount must be non-negative")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = (signal_values + pseudocount) / (background + pseudocount)
        if operation == "ratio":
            return cast(np.ndarray[Any, Any], ratio)
        if operation == "log2_ratio":
            return cast(np.ndarray[Any, Any], np.log2(ratio))
    raise ValueError("operation must be subtract, ratio, or log2_ratio")
