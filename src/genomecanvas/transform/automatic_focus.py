"""Deterministic automatic highlight and zoom selection."""

from __future__ import annotations

import numpy as np

from ..core.intervals import GenomicInterval
from ..model.track import SignalSlice


def largest_peak(signal: SignalSlice, width: int) -> GenomicInterval:
    """Return a fixed-width interval centered on the largest finite value."""
    if width < 1:
        raise ValueError("focus width must be positive")
    values = np.asarray(signal.values, dtype=float)
    if not np.isfinite(values).any():
        raise ValueError("cannot select a focus interval from empty signal")
    index = int(np.nanargmax(values))
    center = int((signal.bin_starts[index] + signal.bin_ends[index]) / 2)
    start = max(signal.start, center - width // 2)
    return GenomicInterval(signal.chrom, start, min(signal.end, start + width))


def highest_enrichment(signal: SignalSlice, background: SignalSlice, width: int) -> GenomicInterval:
    """Select the strongest pointwise signal/background enrichment."""
    if len(signal.values) != len(background.values):
        raise ValueError("signal and background bins must match")
    with np.errstate(divide="ignore", invalid="ignore"):
        enrichment = signal.values / background.values
    selected = SignalSlice(
        signal.chrom, signal.start, signal.end, signal.bin_starts, signal.bin_ends,
        enrichment, signal.units,
    )
    return largest_peak(selected, width)
