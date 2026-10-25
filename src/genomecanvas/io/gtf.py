"""GTF reader with explicit 1-based closed to internal conversion."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from ..core.coordinates import chrom_candidates
from ..core.intervals import GenomicInterval
from ..model.transcript import Transcript


def _attributes(value: str) -> dict[str, str]:
    return {key: val for key, val in re.findall(r'(\w+)\s+"([^"]*)"', value)}


def read_gtf(path: str | Path, *, chrom: str | None = None) -> list[Transcript]:
    """Read transcript/exon records from a GTF file."""
    records: dict[str, dict[str, Any]] = {}
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if len(fields) != 9:
            raise ValueError(f"GTF line {line_number} must contain nine columns")
        record_chrom, _, feature_type, raw_start, raw_end, _, strand, _, raw_attrs = fields
        if chrom is not None:
            if (
                record_chrom not in chrom_candidates(chrom)
                and chrom not in chrom_candidates(record_chrom)
            ):
                continue
        attrs = _attributes(raw_attrs)
        transcript_id = attrs.get("transcript_id")
        gene_id = attrs.get("gene_id", "unknown")
        gene_name = attrs.get("gene_name") or attrs.get("gene") or gene_id
        if transcript_id is None:
            continue
        start, end = int(raw_start) - 1, int(raw_end)
        bucket = records.setdefault(
            transcript_id,
            {"gene_id": gene_id, "gene_name": gene_name, "chrom": record_chrom, "strand": strand,
             "start": start, "end": end, "exons": []},
        )
        bucket["start"] = min(int(bucket["start"]), start)
        bucket["end"] = max(int(bucket["end"]), end)
        if feature_type == "exon":
            bucket["exons"].append(
                GenomicInterval(record_chrom, start, end, strand if strand in {"+", "-"} else None)
            )
    return [
        Transcript(
            transcript_id=transcript_id, gene_id=str(bucket["gene_id"]),
            interval=GenomicInterval(
                str(bucket["chrom"]), int(bucket["start"]), int(bucket["end"]),
                str(bucket["strand"]),
            ),
            exons=tuple(bucket["exons"]),
            name=str(bucket.get("gene_name") or bucket["gene_id"]),
            metadata={"gene_name": str(bucket.get("gene_name") or bucket["gene_id"])},
        )
        for transcript_id, bucket in sorted(records.items())
    ]
