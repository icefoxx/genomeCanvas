"""Optional UCSC Genome Browser REST queries for gene models and cytobands."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Any

from ..core.coordinates import chrom_candidates
from ..core.intervals import GenomicInterval
from ..model.transcript import Transcript
from .cytoband import Cytoband

UCSC_API = "https://api.genome.ucsc.edu"
USER_AGENT = "GenomeCanvas/0.1 (+https://github.com/genomecanvas)"


def _get_json(url: str, timeout: float = 30.0) -> Any:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
        return json.loads(response.read().decode("utf-8"))


def fetch_refgene(
    *,
    build: str,
    chrom: str,
    start: int,
    end: int,
    track: str = "refGene",
) -> list[Transcript]:
    """Download overlapping RefSeq-style transcripts from the UCSC REST API."""
    matched = chrom_candidates(chrom)[0]
    query = urllib.parse.urlencode(
        {
            "genome": build,
            "track": track,
            "chrom": matched,
            "start": str(start),
            "end": str(end),
        }
    )
    payload = _get_json(f"{UCSC_API}/getData/track?{query}")
    rows = payload.get(track, payload.get(matched, []))
    if isinstance(rows, dict):
        rows = rows.get(matched, [])
    transcripts: list[Transcript] = []
    for row in rows:
        name = str(row.get("name") or row.get("transcript") or "unknown")
        gene = str(row.get("name2") or row.get("geneName") or name)
        strand = row.get("strand")
        tx_start = int(row.get("txStart", row.get("start", start)))
        tx_end = int(row.get("txEnd", row.get("end", end)))
        if tx_end <= tx_start:
            continue
        exon_starts = _split_csv(row.get("exonStarts"))
        exon_ends = _split_csv(row.get("exonEnds"))
        exons = tuple(
            GenomicInterval(matched, left, right, strand if strand in {"+", "-"} else None)
            for left, right in zip(exon_starts, exon_ends, strict=False)
            if right > left
        )
        transcripts.append(
            Transcript(
                transcript_id=name,
                gene_id=gene,
                interval=GenomicInterval(
                    matched, tx_start, tx_end, strand if strand in {"+", "-"} else None
                ),
                exons=exons,
                name=gene,
            )
        )
    return transcripts


def fetch_cytobands(*, build: str, chrom: str, track: str = "cytoBandIdeo") -> list[Cytoband]:
    """Download one chromosome of UCSC cytoband ideogram records."""
    matched = chrom_candidates(chrom)[0]
    query = urllib.parse.urlencode({"genome": build, "track": track, "chrom": matched})
    payload = _get_json(f"{UCSC_API}/getData/track?{query}")
    rows = payload.get(track, payload.get(matched, []))
    if isinstance(rows, dict):
        rows = rows.get(matched, [])
    bands: list[Cytoband] = []
    for row in rows:
        start = int(row.get("chromStart", row.get("start", 0)))
        end = int(row.get("chromEnd", row.get("end", 0)))
        if end <= start:
            continue
        bands.append(
            Cytoband(
                str(row.get("chrom", matched)),
                start,
                end,
                str(row.get("name", "")),
                str(row.get("gieStain", row.get("stain", "gneg"))),
            )
        )
    if not bands:
        raise ValueError(f"UCSC returned no cytobands for {build} {chrom}")
    return bands


def fetch_gene_by_symbol(*, build: str, gene: str, track: str = "refGene") -> list[Transcript]:
    """Look up a gene symbol through UCSC refGene and return matching transcripts."""
    query = urllib.parse.urlencode({"genome": build, "track": track})
    payload = _get_json(f"{UCSC_API}/getData/track?{query}")
    rows = payload.get(track, [])
    if isinstance(rows, dict):
        collected: list[Any] = []
        for value in rows.values():
            if isinstance(value, list):
                collected.extend(value)
        rows = collected
    matches: list[Transcript] = []
    needle = gene.casefold()
    for row in rows:
        name2 = str(row.get("name2") or "")
        name = str(row.get("name") or "")
        if needle not in {name2.casefold(), name.casefold()}:
            continue
        chrom = str(row.get("chrom") or "")
        tx_start = int(row.get("txStart", 0))
        tx_end = int(row.get("txEnd", 0))
        if not chrom or tx_end <= tx_start:
            continue
        strand = row.get("strand")
        matches.append(
            Transcript(
                transcript_id=name or gene,
                gene_id=name2 or gene,
                interval=GenomicInterval(
                    chrom, tx_start, tx_end, strand if strand in {"+", "-"} else None
                ),
                name=name2 or gene,
            )
        )
    if not matches:
        raise ValueError(f"UCSC refGene has no transcripts named {gene!r} on {build}")
    return matches


def _split_csv(value: Any) -> list[int]:
    if value is None:
        return []
    text = str(value).strip().strip(",")
    if not text:
        return []
    return [int(item) for item in text.split(",") if item]