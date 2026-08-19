"""Lightweight input inspection command."""

from pathlib import Path


def inspect_file(path: str) -> dict[str, object]:
    """Report format guess and basic file metadata without full indexing."""
    suffix = Path(path).suffix.lower()
    formats = {
        ".bw": "bigwig", ".bedgraph": "bedgraph", ".bed": "bed", ".gtf": "gtf",
        ".gff": "gff", ".vcf": "vcf",
    }
    return {
        "path": path, "format": formats.get(suffix, "unknown"),
        "size": Path(path).stat().st_size,
    }
