"""Summarize profile matrices across regions and optional sample conditions."""

from __future__ import annotations

from typing import Literal

import numpy as np

from .extract import ProfileMatrix

Stat = Literal["mean", "median", "max", "min"]


def profile_summarize(
    matrices: list[ProfileMatrix],
    *,
    stat: Stat = "mean",
    conditions: dict[str, str] | None = None,
) -> dict[str, np.ndarray]:  # type: ignore[type-arg]
    """Reduce each matrix (or condition group) to one vector across bins."""
    grouped: dict[str, list[np.ndarray]] = {}  # type: ignore[type-arg]
    for matrix in matrices:
        label = conditions.get(matrix.sample, matrix.sample) if conditions else matrix.sample
        grouped.setdefault(label, []).append(_reduce_rows(matrix.values, stat))
    return {label: _reduce_stack(rows, stat) for label, rows in grouped.items()}


def _reduce_rows(values: np.ndarray, stat: Stat) -> np.ndarray:  # type: ignore[type-arg]
    if values.size == 0:
        return np.array([], dtype=float)
    return _reduce_stack(list(values), stat)


def _reduce_stack(rows: list[np.ndarray], stat: Stat) -> np.ndarray:  # type: ignore[type-arg]
    stacked = np.vstack(rows) if rows else np.zeros((0, 0))
    if stacked.size == 0:
        return np.array([], dtype=float)
    if stat == "median":
        reduced = np.nanmedian(stacked, axis=0)
    elif stat == "max":
        reduced = np.nanmax(stacked, axis=0)
    elif stat == "min":
        reduced = np.nanmin(stacked, axis=0)
    else:
        reduced = np.nanmean(stacked, axis=0)
    return np.asarray(reduced, dtype=float)
