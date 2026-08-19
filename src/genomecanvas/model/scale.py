"""Explicit signal scale calculations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class ScaleSpec:
    """A named scale policy."""

    mode: str = "per_track"
    ymin: float | None = None
    ymax: float | None = None
    quantile: float = 0.99
    mad_multiplier: float = 3.0

    def resolve(self, values: np.ndarray) -> tuple[float, float]:  # type: ignore[type-arg]
        """Resolve finite values according to the explicit policy."""
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            raise ValueError("cannot resolve a scale from NaN/inf-only signal")
        if self.mode == "fixed":
            if self.ymin is None or self.ymax is None or self.ymax <= self.ymin:
                raise ValueError("fixed scale requires ymax > ymin")
            return self.ymin, self.ymax
        if self.mode in {"robust", "shared_robust"}:
            if not 0 < self.quantile <= 1:
                raise ValueError("robust scale quantile must be in (0, 1]")
            low = float(np.quantile(finite, 1 - self.quantile))
            high = float(np.quantile(finite, self.quantile))
        elif self.mode == "mad":
            median = float(np.median(finite))
            mad = float(np.median(np.abs(finite - median)))
            low = median - self.mad_multiplier * mad
            high = median + self.mad_multiplier * mad
        else:
            low, high = float(np.min(finite)), float(np.max(finite))
        if self.mode == "symmetric":
            bound = max(abs(low), abs(high))
            return -bound, bound
        if high <= low:
            high = low + 1.0
        return low, high
