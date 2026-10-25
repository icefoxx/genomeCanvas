"""BED interval reader."""

from pathlib import Path

from ..core.intervals import GenomicInterval
from ..model.feature import Feature


def read_bed(
    path: str | Path,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
) -> list[Feature]:
    """Read BED, narrowPeak, and broadPeak-like rows."""
    from ..core.coordinates import chrom_candidates

    features: list[Feature] = []
    accepted = set(chrom_candidates(chrom)) if chrom is not None else None
    view = None
    if chrom is not None and start is not None and end is not None:
        view = GenomicInterval(chrom, start, end)
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith(("track", "browser")):
            continue
        fields = raw.split("\t") if "\t" in raw else raw.split()
        if len(fields) < 3:
            raise ValueError(f"BED line {line_number} must have at least three columns")
        strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else None
        interval = GenomicInterval(fields[0], int(fields[1]), int(fields[2]), strand)
        score = float(fields[4]) if len(fields) > 4 and fields[4] not in {".", "-"} else None
        if accepted is not None and interval.chrom not in accepted:
            continue
        if view is not None and not interval.overlaps(view):
            continue
        features.append(
            Feature(interval, "interval", fields[3] if len(fields) > 3 else None, score)
        )
    return features
