"""Motif occurrence track; scanning is intentionally out of scope."""

from ..io.bed import read_bed
from ..model.feature import Feature


class MotifTrack:
    """Read BED-like motif hits."""

    type = "motif"

    def __init__(self, source: str) -> None:
        self.source = source

    def features(self) -> list[Feature]:
        return read_bed(self.source)
