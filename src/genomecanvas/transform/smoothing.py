"""Deterministic moving-average smoothing."""

import numpy as np


def moving_average(values: np.ndarray, window: int) -> np.ndarray:  # type: ignore[type-arg]
    """Smooth finite signal values with a centered rolling window."""
    if window < 1:
        raise ValueError("window must be positive")
    kernel = np.ones(window) / window
    return np.convolve(values, kernel, mode="same")
