"""Minimal VCF record reader with explicit POS conversion."""

from pathlib import Path


def read_vcf(path: str | Path) -> list[dict[str, str]]:
    """Return small variant dictionaries for SNVs and short indels."""
    records = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        chrom, pos, identifier, ref, alt, *_ = raw.split("\t")
        records.append(
            {"chrom": chrom, "start": str(int(pos) - 1), "id": identifier,
             "ref": ref, "alt": alt}
        )
    return records
