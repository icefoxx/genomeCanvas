"""Resolve a gene symbol to a genomic region from local transcript models."""

from __future__ import annotations

from pathlib import Path

from ..core.intervals import GenomicInterval
from ..model.transcript import Transcript
from .gtf import read_gtf


def transcripts_named(transcripts: list[Transcript], gene: str) -> list[Transcript]:
    """Return transcripts whose gene, transcript, or display name matches ``gene``."""
    needle = gene.casefold()
    matches: list[Transcript] = []
    for transcript in transcripts:
        names = {
            transcript.gene_id.casefold(),
            transcript.transcript_id.casefold(),
            (transcript.name or "").casefold(),
            transcript.metadata.get("gene_name", "").casefold(),
        }
        if needle in names:
            matches.append(transcript)
    return matches


def region_for_gene(
    transcripts: list[Transcript],
    gene: str,
    *,
    padding: int = 0,
) -> GenomicInterval:
    """Return the padded union of matching transcripts."""
    matches = transcripts_named(transcripts, gene)
    if not matches:
        available = sorted({item.gene_id for item in transcripts[:20]})
        raise ValueError(f"no transcripts named {gene!r}; examples: {available}")
    chroms = {item.interval.chrom for item in matches}
    if len(chroms) > 1:
        matches = [item for item in matches if item.interval.chrom == matches[0].interval.chrom]
    chrom = matches[0].interval.chrom
    start = max(0, min(item.interval.start for item in matches) - padding)
    end = max(item.interval.end for item in matches) + padding
    if end <= start:
        raise ValueError(f"resolved gene region for {gene!r} is empty")
    return GenomicInterval(chrom, start, end, matches[0].interval.strand)


def region_from_gtf(path: str | Path, gene: str, *, padding: int = 0) -> GenomicInterval:
    """Load a GTF and resolve ``gene`` to a region."""
    return region_for_gene(read_gtf(path), gene, padding=padding)


def collapse_transcripts(transcripts: list[Transcript]) -> list[Transcript]:
    """Merge transcripts that share a gene id into one exon union per gene."""
    groups: dict[str, list[Transcript]] = {}
    for transcript in transcripts:
        groups.setdefault(transcript.gene_id, []).append(transcript)
    collapsed: list[Transcript] = []
    for gene_id, items in groups.items():
        chrom = items[0].interval.chrom
        strand = items[0].interval.strand
        start = min(item.interval.start for item in items)
        end = max(item.interval.end for item in items)
        exons = [exon for item in items for exon in item.exons]
        collapsed.append(
            Transcript(
                transcript_id=gene_id,
                gene_id=gene_id,
                interval=GenomicInterval(chrom, start, end, strand),
                exons=_merge_exons(exons),
                name=items[0].name or gene_id,
            )
        )
    return collapsed


def overlapping_transcripts(
    transcripts: list[Transcript], chrom: str, start: int, end: int
) -> list[Transcript]:
    """Keep transcripts that overlap a view window."""
    view = GenomicInterval(chrom, start, end)
    return [item for item in transcripts if item.interval.overlaps(view)]


def _merge_exons(exons: list[GenomicInterval]) -> tuple[GenomicInterval, ...]:
    if not exons:
        return ()
    ordered = sorted(exons, key=lambda item: (item.chrom, item.start, item.end))
    merged = [ordered[0]]
    for exon in ordered[1:]:
        previous = merged[-1]
        if exon.chrom == previous.chrom and exon.start <= previous.end:
            merged[-1] = GenomicInterval(
                previous.chrom, previous.start, max(previous.end, exon.end), previous.strand
            )
        else:
            merged.append(exon)
    return tuple(merged)
