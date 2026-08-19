"""Deterministic resolution-aware signal aggregation."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any, Literal, cast

import numpy as np

Aggregation = Literal["mean", "max", "min", "median"]


def bin_signal(
    starts: Iterable[int], ends: Iterable[int], values: Iterable[float], *, region_start: int,
    region_end: int, pixels: int, aggregation: Aggregation = "mean",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    """Aggregate piecewise signal into approximately one bin per pixel."""
    if region_end <= region_start or pixels < 1:
        raise ValueError("region_end must exceed region_start and pixels must be positive")
    source_starts = np.asarray(list(starts), dtype=int)
    source_ends = np.asarray(list(ends), dtype=int)
    source_values = np.asarray(list(values), dtype=float)
    if not (len(source_starts) == len(source_ends) == len(source_values)):
        raise ValueError("signal arrays must have equal lengths")
    count = max(1, min(pixels, region_end - region_start))
    edges = np.linspace(region_start, region_end, count + 1, dtype=int)
    out = np.full(count, np.nan)
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        mask = (source_starts < right) & (source_ends > left) & np.isfinite(source_values)
        selected = source_values[mask]
        if selected.size:
            operation = cast(
                Callable[[Any], Any],
                {"mean": np.mean, "max": np.max, "min": np.min, "median": np.median}[aggregation],
            )
            out[index] = float(operation(selected))
    return edges[:-1], edges[1:], out
