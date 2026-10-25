"""Genome-wide bin summaries (multiBigwigSummary bins replacement)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..io.bigwig import Engine, bigwig_chroms, read_bwtool_region_sums, resolve_engine
from ..io.signal import query_signal
from .metadata import SampleMeta


def write_genomic_bins(chrom_sizes: dict[str, int], size: int, path: str | Path) -> Path:
    """Write non-overlapping genomic bins as a BED file."""
    if size < 1:
        raise ValueError("genomic bin size must be positive")
    destination = Path(path)
    lines: list[str] = []
    for chrom, length in chrom_sizes.items():
        start = 0
        index = 0
        while start + size <= length:
            lines.append(f"{chrom}\t{start}\t{start + size}\tbin{index}")
            start += size
            index += 1
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination


def summarize_bins(
    samples: list[SampleMeta],
    bin_bed: str | Path,
    *,
    engine: Engine = "auto",
) -> tuple[tuple[str, ...], np.ndarray]:  # type: ignore[type-arg]
    """Return labels and a bins×samples mean table."""
    resolved = resolve_engine(engine, batch=True)
    labels = tuple(item.label for item in samples)
    columns = []
    for item in samples:
        if resolved == "bwtool" and item.path.lower().endswith((".bw", ".bigwig")):
            columns.append(read_bwtool_region_sums(item.path, bin_bed))
        else:
            from ..io.bed import read_bed

            features = read_bed(bin_bed)
            values = [
                float(np.nanmean(query_signal(
                    item.path, feature.interval.chrom, feature.interval.start, feature.interval.end,
                    pixels=max(1, feature.interval.length), engine=engine,
                ).values))
                for feature in features
            ]
            columns.append(np.asarray(values, dtype=float))
    table = np.column_stack(columns) if columns else np.zeros((0, 0))
    return labels, np.nan_to_num(table, nan=0.0)


def bins_from_first_sample(samples: list[SampleMeta], size: int, path: str | Path) -> Path:
    """Build a bin BED from the first sample's chromosome sizes."""
    sizes = bigwig_chroms(samples[0].path)
    return write_genomic_bins(sizes, size, path)
