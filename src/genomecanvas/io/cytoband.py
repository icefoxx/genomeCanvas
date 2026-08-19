"""Cytoband ideogram reader."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class Cytoband:
    """One cytogenetic band in 0-based half-open coordinates."""

    chrom: str
    start: int
    end: int
    name: str
    stain: str


def read_cytobands(path: str | Path, *, chrom: str) -> list[Cytoband]:
    """Read one chromosome from a UCSC-style cytoband file."""
    result: list[Cytoband] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[0] == chrom:
                result.append(
                    Cytoband(fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4])
                )
    if not result:
        raise ValueError(f"no cytobands found for chromosome {chrom!r}")
    return result
