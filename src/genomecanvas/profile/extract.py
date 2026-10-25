"""Extract per-region signal matrices from BigWig, BedGraph, or peak files."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from ..core.intervals import GenomicInterval
from ..io.bed import read_bed
from ..io.bigwig import Engine, python_profile_matrix, read_bwtool_matrix, resolve_engine
from ..io.genes import overlapping_transcripts
from ..io.gtf import read_gtf
from ..io.signal import infer_signal_format, query_signal
from ..model.transcript import Transcript

StartFrom = Literal["start", "center", "end"]


@dataclass(frozen=True, slots=True)
class ProfileMatrix:
    """One sample's region-by-bin signal matrix."""

    sample: str
    values: np.ndarray  # type: ignore[type-arg]
    up: int
    down: int
    binsize: int
    start_from: StartFrom
    region_ids: tuple[str, ...]
    source: str

    def __post_init__(self) -> None:
        if self.values.ndim != 2:
            raise ValueError("profile matrix must be two-dimensional")
        if self.values.shape[0] != len(self.region_ids):
            raise ValueError("profile matrix rows must match region_ids")


def windows_from_bed(
    path: str | Path,
    *,
    start_from: StartFrom = "center",
    up: int = 2500,
    down: int = 2500,
) -> list[tuple[str, GenomicInterval]]:
    """Build strand-aware windows around BED intervals."""
    windows: list[tuple[str, GenomicInterval]] = []
    for index, feature in enumerate(read_bed(path)):
        interval = _anchor_window(feature.interval, start_from=start_from, up=up, down=down)
        name = feature.name or f"region{index + 1}"
        windows.append((name, interval))
    return windows


def windows_from_transcripts(
    transcripts: list[Transcript],
    *,
    start_from: StartFrom = "start",
    up: int = 2500,
    down: int = 2500,
) -> list[tuple[str, GenomicInterval]]:
    """Build windows around transcript TSS or TES."""
    windows: list[tuple[str, GenomicInterval]] = []
    for transcript in transcripts:
        interval = _anchor_window(transcript.interval, start_from=start_from, up=up, down=down)
        windows.append((transcript.transcript_id, interval))
    return windows


def extract_profile(
    sources: dict[str, str],
    *,
    bed: str | Path | None = None,
    gtf: str | Path | None = None,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    start_from: StartFrom = "center",
    up: int = 2500,
    down: int = 2500,
    binsize: int = 50,
    input_type: Literal["auto", "bw", "bedgraph", "peak"] = "auto",
    engine: Engine = "auto",
) -> list[ProfileMatrix]:
    """Extract a profile matrix for each named signal source."""
    if bed is not None:
        windows = windows_from_bed(bed, start_from=start_from, up=up, down=down)
        region_ids = tuple(name for name, _ in windows)
        bwtool_bed = Path(bed)
    elif gtf is not None:
        transcripts = read_gtf(gtf, chrom=chrom)
        if chrom is not None and start is not None and end is not None:
            transcripts = overlapping_transcripts(transcripts, chrom, start, end)
        windows = windows_from_transcripts(transcripts, start_from=start_from, up=up, down=down)
        region_ids = tuple(name for name, _ in windows)
        bwtool_bed = None
    else:
        raise ValueError("extract_profile requires a BED file or a GTF gene model")
    n_bins = max(1, (up + down) // max(1, binsize))
    coords = [(interval.chrom, interval.start, interval.end) for _, interval in windows]
    resolved = resolve_engine(engine, batch=True)
    matrices: list[ProfileMatrix] = []
    for sample, path in sources.items():
        kind = infer_signal_format(path, input_type)
        if kind == "bw" and resolved == "bwtool":
            if bwtool_bed is not None:
                values = read_bwtool_matrix(
                    path,
                    bwtool_bed,
                    up=up,
                    down=down,
                    binsize=binsize,
                    start_from=start_from,
                )
            else:
                with tempfile.TemporaryDirectory(prefix="gcanvas-anchors-") as tmp:
                    anchor_bed = _write_anchor_bed(windows, Path(tmp))
                    values = read_bwtool_matrix(
                        path,
                        anchor_bed,
                        up=up,
                        down=down,
                        binsize=binsize,
                        start_from="center",
                    )
        elif kind == "bw":
            values = _orient(
                python_profile_matrix(path, coords, n_bins=n_bins, engine=resolved), windows
            )
        else:
            rows = [
                query_signal(
                    path,
                    interval.chrom,
                    interval.start,
                    interval.end,
                    pixels=n_bins,
                    input_type=input_type,
                    engine=engine,
                ).values
                for _, interval in windows
            ]
            values = _orient(np.vstack(rows) if rows else np.zeros((0, n_bins)), windows)
        matrices.append(
            ProfileMatrix(
                sample=sample,
                values=values,
                up=up,
                down=down,
                binsize=binsize,
                start_from=start_from,
                region_ids=region_ids,
                source=path,
            )
        )
    return matrices


def _orient(values: np.ndarray, windows: list[tuple[str, GenomicInterval]]) -> np.ndarray:  # type: ignore[type-arg]
    """Reverse minus-strand rows so upstream is always on the left (as deepTools and bwtool do)."""
    minus = np.array([interval.strand == "-" for _, interval in windows], dtype=bool)
    if values.size and minus.any():
        values = np.array(values, copy=True)
        values[minus] = values[minus, ::-1]
    return values


def _write_anchor_bed(windows: list[tuple[str, GenomicInterval]], directory: Path) -> Path:
    """Write 1 bp anchors so bwtool matrix expands from the same points."""
    path = Path(directory) / "gcanvas-anchors.bed"
    lines = []
    for name, interval in windows:
        # windows are already expanded; the midpoint is the original anchor.
        anchor = (interval.start + interval.end) // 2
        strand = interval.strand or "."
        lines.append(f"{interval.chrom}\t{anchor}\t{anchor + 1}\t{name}\t0\t{strand}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _anchor_window(
    interval: GenomicInterval, *, start_from: StartFrom, up: int, down: int
) -> GenomicInterval:
    strand = interval.strand
    if start_from == "center":
        # Round toward the 5' end on each strand so +/- windows are exact mirror images.
        anchor = (interval.start + interval.end + (strand == "-")) // 2
    elif start_from == "end":
        anchor = interval.end if strand != "-" else interval.start
    else:
        anchor = interval.start if strand != "-" else interval.end
    left, right = up, down
    if strand == "-":
        left, right = down, up
    start = max(0, anchor - left)
    end = max(start + 1, anchor + right)
    return GenomicInterval(interval.chrom, start, end, strand)
