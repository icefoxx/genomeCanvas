"""Immutable 0-based, half-open genomic intervals."""

from __future__ import annotations

from dataclasses import dataclass, replace

from .errors import CoordinateError


@dataclass(frozen=True, slots=True)
class GenomicInterval:
    """A validated interval using the internal convention ``[start, end)``."""

    chrom: str
    start: int
    end: int
    strand: str | None = None

    def __post_init__(self) -> None:
        if not self.chrom:
            raise CoordinateError("chrom must not be empty")
        if isinstance(self.start, bool) or isinstance(self.end, bool):
            raise CoordinateError("start and end must be integers")
        if self.start < 0:
            raise CoordinateError("start must be >= 0")
        if self.end <= self.start:
            raise CoordinateError("end must be greater than start")
        if self.strand not in (None, "+", "-"):
            raise CoordinateError("strand must be '+', '-', or None")

    @property
    def length(self) -> int:
        """Return interval width in base pairs."""
        return self.end - self.start

    def overlaps(self, other: GenomicInterval) -> bool:
        """Return whether intervals overlap on the same chromosome."""
        return self.chrom == other.chrom and self.start < other.end and other.start < self.end

    def intersection(self, other: GenomicInterval) -> GenomicInterval | None:
        """Return the overlap, or ``None`` when there is no overlap."""
        if not self.overlaps(other):
            return None
        return GenomicInterval(self.chrom, max(self.start, other.start), min(self.end, other.end))

    def contains(self, position: int | GenomicInterval) -> bool:
        """Return whether a base position or complete interval is contained."""
        if isinstance(position, GenomicInterval):
            return (
                self.chrom == position.chrom
                and self.start <= position.start
                and position.end <= self.end
            )
        return self.start <= position < self.end

    def distance_to(self, other: GenomicInterval) -> int:
        """Return gap distance, or zero for overlapping intervals."""
        if self.chrom != other.chrom:
            raise CoordinateError("distance is undefined for different chromosomes")
        if self.overlaps(other):
            return 0
        return max(other.start - self.end, self.start - other.end)

    def expand(self, left: int = 0, right: int | None = None) -> GenomicInterval:
        """Expand by base pairs without allowing the start below zero."""
        if left < 0 or (right is not None and right < 0):
            raise CoordinateError("expansion widths must be non-negative")
        right = left if right is None else right
        return replace(self, start=max(0, self.start - left), end=self.end + right)

    def clip(self, start: int, end: int) -> GenomicInterval | None:
        """Clip to a valid half-open view, preserving original interval semantics."""
        if start < 0 or end <= start:
            raise CoordinateError("clip bounds must satisfy 0 <= start < end")
        clipped_start = max(self.start, start)
        clipped_end = min(self.end, end)
        if clipped_end <= clipped_start:
            return None
        return replace(self, start=clipped_start, end=clipped_end)
