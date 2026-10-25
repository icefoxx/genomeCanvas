"""BigWig metadata TSV used by the heatmap pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

COLUMNS = ("bigWig", "bigWigLabel", "conditionGrp", "metaGrp", "colorMap", "profileColor")


@dataclass(frozen=True, slots=True)
class SampleMeta:
    """One BigWig sample and its plot annotations."""

    path: str
    label: str
    condition: str
    meta_group: str
    colormap: str
    profile_color: str


def read_sample_metadata(
    path: str | Path,
    *,
    default_colormap: str = "Spectral_r",
    default_profile_color: str = "auto",
) -> list[SampleMeta]:
    """Read the six-column heatmap metadata TSV (header optional)."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"metadata file does not exist: {source}")
    rows: list[list[str]] = []
    saw_header = False
    for raw in source.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        fields = [item.strip() for item in raw.rstrip("\n\r").split("\t")]
        if _is_header(fields):
            if saw_header:
                continue
            saw_header = True
            continue
        if line.lstrip().startswith("#"):
            continue
        rows.append(fields)
    if not rows:
        raise ValueError("metadata file has no data rows")
    samples: list[SampleMeta] = []
    labels: list[str] = []
    for fields in rows:
        padded = list(fields)
        while len(padded) < 6:
            padded.append("")
        if len(padded) > 6:
            padded = padded[:6]
        if not padded[2]:
            padded[2] = "def"
        if not padded[3]:
            padded[3] = "def"
        if not padded[4]:
            padded[4] = default_colormap
        if not padded[5]:
            padded[5] = default_profile_color
        label = padded[1] or Path(padded[0]).stem
        if label in labels:
            raise ValueError(f"duplicate sample label {label!r}")
        labels.append(label)
        bw_path = padded[0]
        resolved = Path(bw_path)
        if not resolved.is_file():
            resolved = (source.parent / bw_path).resolve()
            bw_path = str(resolved)
        samples.append(
            SampleMeta(
                path=bw_path,
                label=label,
                condition=padded[2],
                meta_group=padded[3],
                colormap=padded[4],
                profile_color=padded[5],
            )
        )
    missing = [item.path for item in samples if not Path(item.path).is_file()]
    if missing:
        raise FileNotFoundError("BigWig files are not accessible: " + ", ".join(missing))
    return samples


def _is_header(fields: list[str]) -> bool:
    if len(fields) != 6:
        return False
    first = fields[0].lstrip("#")
    return [first, *fields[1:]] == list(COLUMNS)
