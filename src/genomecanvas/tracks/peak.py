"""Peak interval track."""

from ..io.bed import read_bed
from ..model.feature import Feature


class PeakTrack:
    """Read BED/narrowPeak intervals without interpreting caller scores."""

    type = "peak"

    def __init__(self, source: str) -> None:
        self.source = source

    def features(self) -> list[Feature]:
        """Return source features."""
        return read_bed(self.source)
