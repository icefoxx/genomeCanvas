"""Unified regional signal queries for BedGraph, BigWig, and peak files."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np

from ..model.track import SignalSlice
from ..transform.binning import Aggregation, bin_signal
from .bed import read_bed
from .bedgraph import read_bedgraph
from .bigwig import Engine, binned_bigwig

SignalFormat = Literal["auto", "bw", "bedgraph", "peak"]


def infer_signal_format(path: str | Path, input_type: SignalFormat = "auto") -> SignalFormat:
    """Infer the on-disk signal format from an explicit type or file suffix."""
    if input_type != "auto":
        return input_type
    suffix = Path(path).suffix.lower()
    if suffix in {".bw", ".bigwig", ".bigwig.gz"}:
        return "bw"
    if suffix in {".narrowpeak", ".broadpeak"}:
        return "peak"
    name = Path(path).name.lower()
    if name.endswith(".narrowpeak.gz") or name.endswith(".broadpeak.gz") or "narrowpeak" in name:
        return "peak"
    return "bedgraph"


def query_signal(
    path: str | Path,
    chrom: str,
    start: int,
    end: int,
    *,
    pixels: int = 1000,
    aggregation: Aggregation = "mean",
    input_type: SignalFormat = "auto",
    engine: Engine = "auto",
) -> SignalSlice:
    """Return a binned signal slice for one genomic window."""
    kind = infer_signal_format(path, input_type)
    n_bins = max(1, min(pixels, end - start))
    if kind == "bw":
        bw_stat: Literal["mean", "max", "min", "sum"] = "mean"
        if aggregation in ("max", "min"):
            bw_stat = aggregation
        values = binned_bigwig(
            str(path), chrom, start, end, n_bins=n_bins, aggregation=bw_stat, engine=engine
        )
        edges = np.linspace(start, end, n_bins + 1, dtype=int)
        return SignalSlice(chrom, start, end, edges[:-1], edges[1:], values)
    if kind == "peak":
        features = [
            feature
            for feature in read_bed(path)
            if feature.interval.chrom == chrom
            and feature.interval.start < end
            and feature.interval.end > start
        ]
        starts = [feature.interval.start for feature in features]
        ends = [feature.interval.end for feature in features]
        scores = [float(feature.score or 0.0) for feature in features]
        bin_starts, bin_ends, binned = bin_signal(
            starts, ends, scores, region_start=start, region_end=end, pixels=n_bins,
            aggregation="max",
        )
        return SignalSlice(chrom, start, end, bin_starts, bin_ends, np.nan_to_num(binned, nan=0.0))
    records = read_bedgraph(path, chrom=chrom, start=start, end=end)
    starts = [record.start for record, _ in records]
    ends = [record.end for record, _ in records]
    scores = [value for _, value in records]
    bin_starts, bin_ends, binned = bin_signal(
        starts, ends, scores, region_start=start, region_end=end, pixels=n_bins,
        aggregation=aggregation,
    )
    return SignalSlice(
        chrom, start, end, bin_starts, bin_ends, np.nan_to_num(binned, nan=0.0)
    )
