"""Generic matrix track."""

import numpy as np


class HeatmapTrack:
    """Hold a deterministic dense matrix and genomic bounds."""

    type = "heatmap"

    def __init__(self, matrix: np.ndarray, start: int, end: int) -> None:  # type: ignore[type-arg]
        if matrix.ndim != 2:
            raise ValueError("heatmap matrix must be two-dimensional")
        self.matrix, self.start, self.end = matrix, start, end
