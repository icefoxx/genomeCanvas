"""Signal transformation helpers."""

from typing import cast

import numpy as np


def transform_signal(values: np.ndarray, operation: str) -> np.ndarray:  # type: ignore[type-arg]
    """Apply a named display-only transform."""
    if operation == "identity":
        return values.copy()
    if operation == "log1p":
        return cast(np.ndarray, np.log1p(values))  # type: ignore[type-arg]
    raise ValueError("unsupported signal transform")
