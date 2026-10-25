"""Chromatin-state (chromHMM) BED reader and UCSC Broad HMM colors."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..core.intervals import GenomicInterval
from .bed import read_bed

# Colors follow the public UCSC wgEncodeBroadHmm 15-state legend.
ENCODE_BROAD_HMM_COLORS: dict[str, str] = {
    "1": "#FF0000",
    "2": "#FF4500",
    "3": "#FF4500",
    "4": "#FFC34D",
    "5": "#FFC34D",
    "6": "#FFFF00",
    "7": "#FFFF00",
    "8": "#008000",
    "9": "#006400",
    "10": "#00C800",
    "11": "#3CB371",
    "12": "#7E26A5",
    "13": "#808080",
    "14": "#C0C0C0",
    "15": "#C0C0C0",
}

# UCSC BED score shading used for narrowPeak/broadPeak intensity tracks.
UCSC_SCORE_COLORS: tuple[tuple[float, str], ...] = (
    (166, "#FFFFFF"),
    (277, "#F0F0F0"),
    (389, "#D9D9D9"),
    (500, "#BDBDBD"),
    (612, "#969696"),
    (723, "#737373"),
    (834, "#525252"),
    (945, "#252525"),
    (float("inf"), "#000000"),
)


@dataclass(frozen=True, slots=True)
class ChromatinState:
    """One chromHMM interval with a state label."""

    interval: GenomicInterval
    state: str


def read_chromhmm(path: str | Path) -> list[ChromatinState]:
    """Read a 4-column chromHMM BED (chrom, start, end, state)."""
    states: list[ChromatinState] = []
    for feature in read_bed(path):
        states.append(ChromatinState(feature.interval, feature.name or "1"))
    return states


def chromhmm_color(state: str, mapping: dict[str, str] | None = None) -> str:
    """Resolve a state name or integer code to a display color."""
    colors: dict[str, str] = mapping or ENCODE_BROAD_HMM_COLORS
    key = str(state)
    if key in colors:
        return colors[key]
    digits = "".join(character for character in state if character.isdigit())
    if digits in colors:
        return colors[digits]
    return "#888888"


def ucsc_score_color(score: float) -> str:
    """Map a 0–1000-style score onto the UCSC greyscale convention."""
    for threshold, color in UCSC_SCORE_COLORS:
        if score <= threshold:
            return color
    return "#000000"
