"""Gene model track backed by GTF transcripts."""

from ..io.gtf import read_gtf
from ..model.transcript import Transcript


class GeneTrack:
    """Expose transcript structures for a selected chromosome."""

    type = "gene"

    def __init__(self, source: str) -> None:
        self.source = source

    def transcripts(self, chrom: str | None = None) -> list[Transcript]:
        return read_gtf(self.source, chrom=chrom)
