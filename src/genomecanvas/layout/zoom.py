"""Coordinate mapping for parent/child zoom panels."""

from dataclasses import dataclass

from ..core.intervals import GenomicInterval


@dataclass(frozen=True, slots=True)
class ZoomSpec:
    parent: GenomicInterval
    child: GenomicInterval

    def __post_init__(self) -> None:
        if self.parent.chrom != self.child.chrom or not self.parent.overlaps(self.child):
            raise ValueError("child zoom must overlap parent region")

    def map_position(self, position: int, width: float) -> float:
        """Map a genomic position to child-panel pixels."""
        return (position - self.child.start) * width / self.child.length
