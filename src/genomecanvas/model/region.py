"""Immutable locus-level region model."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from ..core.intervals import GenomicInterval


@dataclass(frozen=True, slots=True)
class GenomeRegion:
    """A validated genomic view with assembly and presentation metadata."""

    chrom: str
    start: int
    end: int
    assembly: str
    label: str | None = None
    gene: str | None = None
    strand: str | None = None
    metadata: tuple[tuple[str, Any], ...] = ()

    def __post_init__(self) -> None:
        GenomicInterval(self.chrom, self.start, self.end, self.strand)
        if not self.assembly:
            raise ValueError("assembly must not be empty")

    @property
    def interval(self) -> GenomicInterval:
        """Return the region as a core interval."""
        return GenomicInterval(self.chrom, self.start, self.end, self.strand)

    @classmethod
    def parse(cls, value: str, *, assembly: str, label: str | None = None) -> GenomeRegion:
        """Parse ``chrom:start-end`` with optional thousands separators."""
        match = re.fullmatch(r"([^:]+):([0-9][0-9,]*)-([0-9][0-9,]*)", value.strip())
        if match is None:
            raise ValueError("region must have the form chrom:start-end")
        chrom, start, end = match.groups()
        return cls(chrom, int(start.replace(",", "")), int(end.replace(",", "")), assembly, label)
