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


def _open_text(path: Path):  # type: ignore[no-untyped-def]
    if path.suffix == ".gz" or path.name.endswith(".txt.gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def stain_color(stain: str) -> str:
    """Map a UCSC gieStain token to an ideogram fill color."""
    mapping = {
        "gneg": "#FFFFFF",
        "gpos25": "#C8C8C8",
        "gpos50": "#969696",
        "gpos75": "#646464",
        "gpos100": "#000000",
        "acen": "#D62728",
        "gvar": "#B0B0B0",
        "stalk": "#6BAED6",
    }
    return mapping.get(stain, "#D0D0D0")


def read_cytobands(path: str | Path, *, chrom: str) -> list[Cytoband]:
    """Read one chromosome from a UCSC-style cytoband file."""
    from ..core.coordinates import chrom_candidates

    result: list[Cytoband] = []
    accepted = set(chrom_candidates(chrom))
    with _open_text(Path(path)) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[0] in accepted:
                result.append(
                    Cytoband(fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4])
                )
    if not result:
        raise ValueError(f"no cytobands found for chromosome {chrom!r}")
    return result
