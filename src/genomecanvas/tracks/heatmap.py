"""Generic matrix track and profile-heatmap adapter."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..profile.extract import ProfileMatrix


class HeatmapTrack:
    """Hold a deterministic dense matrix and genomic bounds."""

    type = "heatmap"

    def __init__(self, matrix: np.ndarray, start: int, end: int) -> None:  # type: ignore[type-arg]
        if matrix.ndim != 2:
            raise ValueError("heatmap matrix must be two-dimensional")
        self.matrix, self.start, self.end = matrix, start, end

    @classmethod
    def from_profile(cls, profile: ProfileMatrix) -> HeatmapTrack:
        """Wrap a profile matrix as a heatmap spanning ``[-up, down]``."""
        return cls(profile.values, -profile.up, profile.down)
