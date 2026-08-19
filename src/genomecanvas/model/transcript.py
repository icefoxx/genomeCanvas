"""Transcript and exon models."""

from dataclasses import dataclass, field

from ..core.intervals import GenomicInterval


@dataclass(frozen=True, slots=True)
class Transcript:
    """A transcript with ordered exon intervals."""

    transcript_id: str
    gene_id: str
    interval: GenomicInterval
    exons: tuple[GenomicInterval, ...] = ()
    cds: tuple[GenomicInterval, ...] = ()
    name: str | None = None
    metadata: dict[str, str] = field(default_factory=dict)
