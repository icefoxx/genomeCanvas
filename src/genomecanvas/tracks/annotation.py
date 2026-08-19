"""Generic interval annotation track."""

from ..io.bed import read_bed
from ..model.feature import Feature


class AnnotationTrack:
    """Display user-supplied interval annotations."""

    type = "annotation"

    def __init__(self, source: str) -> None:
        self.source = source

    def features(self) -> list[Feature]:
        return read_bed(self.source)


class HighlightTrack:
    """Pure visual overlays that do not modify underlying features."""

    type = "highlight"

    def __init__(self, intervals: list[tuple[int, int, str | None]]) -> None:
        self.intervals = tuple(intervals)
