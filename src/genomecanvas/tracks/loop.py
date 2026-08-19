"""Paired-interval loop data and deterministic arc heights."""

from pathlib import Path


class LoopTrack:
    """Read BEDPE-like paired anchors."""

    type = "loop"

    def __init__(self, source: str) -> None:
        self.source = source

    def loops(self) -> list[dict[str, int | str]]:
        result: list[dict[str, int | str]] = []
        for raw in Path(self.source).read_text(encoding="utf-8").splitlines():
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split("\t")
            if len(fields) < 6:
                raise ValueError("BEDPE rows require six columns")
            result.append(
                {"chrom1": fields[0], "start1": int(fields[1]), "end1": int(fields[2]),
                 "chrom2": fields[3], "start2": int(fields[4]), "end2": int(fields[5])}
            )
        return result
