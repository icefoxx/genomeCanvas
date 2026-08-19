"""Explicit visual comparison of two compatible signal slices."""

from __future__ import annotations

import warnings

import numpy as np

from ..model.track import SignalSlice


def compare_signal(
    left: SignalSlice,
    right: SignalSlice,
    *,
    operation: str = "difference",
    pseudocount: float | None = None,
) -> SignalSlice:
    """Compute difference, ratio, or log2 ratio with explicit safeguards."""
    if (left.chrom, left.start, left.end) != (right.chrom, right.start, right.end):
        raise ValueError("signal comparison requires matching chromosome and region")
    if left.units != right.units:
        raise ValueError("signal comparison requires matching declared units")
    if len(left.values) != len(right.values) or not np.array_equal(
        left.bin_starts, right.bin_starts
    ):
        raise ValueError("signal comparison requires matching bins")
    if operation not in {"difference", "ratio", "log2_ratio"}:
        raise ValueError("operation must be difference, ratio, or log2_ratio")
    if operation in {"ratio", "log2_ratio"} and pseudocount is None:
        raise ValueError("ratio operations require an explicit pseudocount")
    if pseudocount is not None and pseudocount < 0:
        raise ValueError("pseudocount must be non-negative")
    if operation == "difference":
        values = right.values - left.values
    else:
        assert pseudocount is not None
        denominator = left.values + pseudocount
        if np.any(denominator == 0):
            warnings.warn(
                "zero denominator produces NaN in signal ratio", RuntimeWarning, stacklevel=2
            )
        with np.errstate(divide="ignore", invalid="ignore"):
            values = right.values + pseudocount
            values = values / denominator
            if operation == "log2_ratio":
                values = np.log2(values)
    return SignalSlice(
        left.chrom, left.start, left.end, left.bin_starts, left.bin_ends, values, left.units
    )
