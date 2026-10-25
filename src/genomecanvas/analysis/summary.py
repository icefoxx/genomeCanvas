"""Per-region area summaries used by PCA and differential peak tests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..io.bed import read_bed
from ..io.bigwig import Engine, read_bwtool_region_sums, resolve_engine
from ..io.signal import infer_signal_format, query_signal


@dataclass(frozen=True, slots=True)
class RegionSummary:
    """Sample-by-region intensity table."""

    regions: tuple[str, ...]
    samples: tuple[str, ...]
    chroms: tuple[str, ...]
    starts: tuple[int, ...]
    ends: tuple[int, ...]
    values: np.ndarray  # type: ignore[type-arg]


def extract_summary(
    sources: dict[str, str],
    bed: str | Path,
    *,
    input_type: str = "auto",
    engine: Engine = "auto",
) -> RegionSummary:
    """Sum signal over each BED interval for every sample."""
    features = read_bed(bed)
    region_ids = tuple(feature.name or f"r{index}" for index, feature in enumerate(features))
    samples = tuple(sources)
    matrix = np.zeros((len(features), len(samples)), dtype=float)
    resolved = resolve_engine(engine, batch=True)
    for column, sample in enumerate(samples):
        path = sources[sample]
        kind = infer_signal_format(path, input_type)  # type: ignore[arg-type]
        if kind == "bw" and resolved == "bwtool":
            sums = read_bwtool_region_sums(path, bed)
            if sums.size != len(features):
                raise ValueError("bwtool summary rows do not match the BED file")
            matrix[:, column] = sums
            continue
        for row, feature in enumerate(features):
            slice_ = query_signal(
                path,
                feature.interval.chrom,
                feature.interval.start,
                feature.interval.end,
                pixels=max(1, feature.interval.length),
                input_type=input_type,  # type: ignore[arg-type]
                engine=engine,
            )
            matrix[row, column] = float(np.nansum(slice_.values))
    return RegionSummary(
        regions=region_ids,
        samples=samples,
        chroms=tuple(feature.interval.chrom for feature in features),
        starts=tuple(feature.interval.start for feature in features),
        ends=tuple(feature.interval.end for feature in features),
        values=matrix,
    )
