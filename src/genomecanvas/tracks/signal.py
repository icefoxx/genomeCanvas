"""Continuous signal track backed by normalized slices."""

from __future__ import annotations

from typing import Any

from ..io.bedgraph import read_bedgraph
from ..model.track import SignalSlice, Track
from ..transform.binning import Aggregation, bin_signal


class SignalTrack(Track):
    """A signal track that queries only the requested region."""

    def __init__(
        self, id: str, label: str, source: str, aggregation: Aggregation = "mean", **kwargs: Any
    ) -> None:
        super().__init__(id=id, type="signal", label=label, source=source, **kwargs)
        self.aggregation = aggregation

    def query(self, chrom: str, start: int, end: int, *, pixels: int = 1000) -> SignalSlice:
        """Read and bin this track's BedGraph source for a region."""
        if self.source is None:
            raise ValueError(f"signal track {self.id!r} has no source")
        records = read_bedgraph(self.source, chrom=chrom, start=start, end=end)
        starts = [record.start for record, _ in records]
        ends = [record.end for record, _ in records]
        values = [value for _, value in records]
        bin_starts, bin_ends, binned = bin_signal(
            starts, ends, values, region_start=start, region_end=end, pixels=pixels,
            aggregation=self.aggregation,
        )
        return SignalSlice(chrom, start, end, bin_starts, bin_ends, binned, self.units)
