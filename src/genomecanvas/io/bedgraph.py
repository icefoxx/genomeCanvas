"""Small, deterministic BedGraph reader."""

from __future__ import annotations

from pathlib import Path

from ..core.intervals import GenomicInterval


def read_bedgraph(
    path: str | Path, *, chrom: str, start: int, end: int
) -> list[tuple[GenomicInterval, float]]:
    """Read overlapping BedGraph records for one region.

    BedGraph is already 0-based half-open; malformed records fail loudly.
    """
    result: list[tuple[GenomicInterval, float]] = []
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 4:
            raise ValueError(f"BedGraph line {line_number} has fewer than four columns")
        record_chrom, raw_start, raw_end, raw_value = fields[:4]
        record = GenomicInterval(record_chrom, int(raw_start), int(raw_end))
        if record.chrom == chrom and record.overlaps(GenomicInterval(chrom, start, end)):
            result.append((record, float(raw_value)))
    return result
