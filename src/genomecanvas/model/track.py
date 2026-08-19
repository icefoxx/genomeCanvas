"""Normalized track data and shared metadata."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True, slots=True)
class SignalSlice:
    """Regional binned signal independent of its source file format."""

    chrom: str
    start: int
    end: int
    bin_starts: np.ndarray  # type: ignore[type-arg]
    bin_ends: np.ndarray  # type: ignore[type-arg]
    values: np.ndarray  # type: ignore[type-arg]
    units: str | None = None
    metadata: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not (len(self.bin_starts) == len(self.bin_ends) == len(self.values)):
            raise ValueError("signal slice arrays must have equal lengths")
        if self.end <= self.start:
            raise ValueError("signal slice must have positive width")


@dataclass(frozen=True, slots=True)
class Track:
    """Common scientific and rendering metadata for a track."""

    id: str
    type: str
    label: str
    source: str | None = None
    height: float = 1.0
    visibility: bool = True
    metadata: dict[str, object] = field(default_factory=dict)
    assay: str | None = None
    condition: str | None = None
    replicate: str | None = None
    sample: str | None = None
    group: str | None = None
    units: str | None = None
    genome_assembly: str | None = None
    scale_group: str | None = None
    z_order: int = 0
    style: dict[str, object] = field(default_factory=dict)
