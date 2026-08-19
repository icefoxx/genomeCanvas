"""GFF3 adapter with explicit 1-based closed conversion."""

from pathlib import Path

from ..core.intervals import GenomicInterval


def read_gff(path: str | Path) -> list[GenomicInterval]:
    """Read the genomic spans from GFF3 records."""
    result = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        strand = fields[6] if fields[6] in {"+", "-"} else None
        result.append(GenomicInterval(fields[0], int(fields[3]) - 1, int(fields[4]), strand))
    return result
