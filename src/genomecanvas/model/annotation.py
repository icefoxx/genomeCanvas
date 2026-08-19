"""Evidence annotation model."""

from dataclasses import dataclass, field

from ..core.intervals import GenomicInterval


@dataclass(frozen=True, slots=True)
class EvidenceItem:
    id: str
    region: GenomicInterval
    evidence_type: str
    source_track_ids: tuple[str, ...] = ()
    label: str | None = None
    metadata: dict[str, object] = field(default_factory=dict)
