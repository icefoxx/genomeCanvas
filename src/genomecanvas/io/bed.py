"""BED interval reader."""

from pathlib import Path

from ..core.intervals import GenomicInterval
from ..model.feature import Feature


def read_bed(path: str | Path) -> list[Feature]:
    """Read BED, narrowPeak, and broadPeak-like rows."""
    features: list[Feature] = []
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if len(fields) < 3:
            raise ValueError(f"BED line {line_number} must have at least three columns")
        interval = GenomicInterval(fields[0], int(fields[1]), int(fields[2]))
        score = float(fields[4]) if len(fields) > 4 and fields[4] not in {".", "-"} else None
        features.append(
            Feature(interval, "interval", fields[3] if len(fields) > 3 else None, score)
        )
    return features
