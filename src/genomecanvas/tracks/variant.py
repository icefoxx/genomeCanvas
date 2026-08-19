"""Minimal VCF variant track."""

from pathlib import Path


class VariantTrack:
    """Read SNVs and short indels without clinical interpretation."""

    type = "variant"

    def __init__(self, source: str) -> None:
        self.source = source

    def variants(self) -> list[dict[str, str]]:
        result = []
        for raw in Path(self.source).read_text(encoding="utf-8").splitlines():
            if not raw.strip() or raw.startswith("#"):
                continue
            chrom, pos, identifier, ref, alt, *_ = raw.split("\t")
            result.append(
                {"chrom": chrom, "pos": str(int(pos) - 1), "id": identifier,
                 "ref": ref, "alt": alt}
            )
        return result
