"""Explicit, user-requested normalization operations."""

import numpy as np


def zscore(values: np.ndarray) -> np.ndarray:  # type: ignore[type-arg]
    """Standardize values; does not infer biological normalization."""
    mean, std = float(np.nanmean(values)), float(np.nanstd(values))
    if std == 0:
        raise ValueError("cannot z-score constant signal")
    return (values - mean) / std
