"""Generic genomic feature records."""

from dataclasses import dataclass, field

from ..core.intervals import GenomicInterval


@dataclass(frozen=True, slots=True)
class Feature:
    """An interval with a typed feature kind and metadata."""

    interval: GenomicInterval
    kind: str
    name: str | None = None
    score: float | None = None
    metadata: dict[str, str] = field(default_factory=dict)
