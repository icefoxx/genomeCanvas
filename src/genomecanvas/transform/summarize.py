"""Collapse aligned signal slices by sample condition."""

from __future__ import annotations

from typing import Literal

import numpy as np

from ..model.track import SignalSlice

Stat = Literal["mean", "median", "max", "min"]


def summarize_signals(
    slices: dict[str, SignalSlice],
    conditions: dict[str, str],
    *,
    stat: Stat = "mean",
) -> dict[str, SignalSlice]:
    """Average (or otherwise reduce) slices that share a condition label."""
    groups: dict[str, list[tuple[str, SignalSlice]]] = {}
    for sample, slice_ in slices.items():
        groups.setdefault(conditions[sample], []).append((sample, slice_))
    result: dict[str, SignalSlice] = {}
    for condition, items in groups.items():
        reference = items[0][1]
        stacked = np.vstack([item.values for _, item in items])
        if stat == "median":
            reduced = np.nanmedian(stacked, axis=0)
        elif stat == "max":
            reduced = np.nanmax(stacked, axis=0)
        elif stat == "min":
            reduced = np.nanmin(stacked, axis=0)
        else:
            reduced = np.nanmean(stacked, axis=0)
        values = np.asarray(reduced, dtype=float)
        result[condition] = SignalSlice(
            reference.chrom,
            reference.start,
            reference.end,
            reference.bin_starts,
            reference.bin_ends,
            values,
            reference.units,
            metadata={"stat": stat, "members": ",".join(name for name, _ in items)},
        )
    return result
