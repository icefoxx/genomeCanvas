"""Amplitude and shape metrics on physical bin coordinates (bp relative to the anchor).

Mass inside a bin is treated as uniformly spread over that bin, so cumulative distributions are
piecewise linear between bin edges. That makes W1, quantiles and asymmetry exact for the
binned representation and valid for non-uniform bin widths.
"""

from __future__ import annotations

import numpy as np

EPS = 1e-12


# ---------------------------------------------------------------- binned distributions


def widths(edges: np.ndarray) -> np.ndarray:
    return np.diff(np.asarray(edges, dtype=float))


def centers(edges: np.ndarray) -> np.ndarray:
    edges = np.asarray(edges, dtype=float)
    return (edges[:-1] + edges[1:]) / 2


def normalise(mass: np.ndarray) -> np.ndarray:
    """Bin masses → probabilities along the last axis (NaN where total mass is not positive)."""
    mass = np.asarray(mass, dtype=float)
    total = mass.sum(axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(total > EPS, mass / np.where(total > EPS, total, 1.0), np.nan)


def cdf_at_edges(p: np.ndarray) -> np.ndarray:
    """F evaluated at every bin edge (prepends 0)."""
    p = np.asarray(p, dtype=float)
    zeros = np.zeros(p.shape[:-1] + (1,))
    return np.concatenate([zeros, np.cumsum(p, axis=-1)], axis=-1)


def centroid(p: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return np.sum(p * centers(edges), axis=-1)


def spread_sd(p: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Standard deviation (bp), including the within-bin variance w²/12 of uniform mass."""
    c, w = centers(edges), widths(edges)
    mu = np.sum(p * c, axis=-1, keepdims=True)
    var = np.sum(p * ((c - mu) ** 2 + w**2 / 12.0), axis=-1)
    return np.sqrt(np.maximum(var, 0.0))


def mean_abs_deviation(p: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Mean absolute deviation from the centroid (bp); less sensitive to window-edge noise than sd."""
    c = centers(edges)
    mu = np.sum(p * c, axis=-1, keepdims=True)
    return np.sum(p * np.abs(c - mu), axis=-1)


def quantile(p: np.ndarray, edges: np.ndarray, q: float) -> np.ndarray:
    """Inverse of the piecewise-linear CDF (first crossing; flat stretches are skipped)."""
    edges = np.asarray(edges, dtype=float)
    shape = np.shape(p)[:-1]
    rows = np.asarray(p, dtype=float).reshape(-1, np.shape(p)[-1])
    F = cdf_at_edges(rows)
    k = np.argmax(F >= q - 1e-12, axis=1)  # first edge where F reaches q
    k = np.clip(k, 1, len(edges) - 1)
    lo = np.take_along_axis(F, (k - 1)[:, None], axis=1)[:, 0]
    hi = np.take_along_axis(F, k[:, None], axis=1)[:, 0]
    frac = np.where(hi - lo > EPS, (q - lo) / np.maximum(hi - lo, EPS), 0.0)
    out = edges[k - 1] + np.clip(frac, 0, 1) * (edges[k] - edges[k - 1])
    out[~np.all(np.isfinite(rows), axis=1)] = np.nan
    return out.reshape(shape)


def iqr_width(p: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return quantile(p, edges, 0.75) - quantile(p, edges, 0.25)


def asymmetry(p: np.ndarray, edges: np.ndarray, origin: float = 0.0) -> np.ndarray:
    """(mass right of origin − mass left of origin) / total, splitting a straddling bin."""
    edges = np.asarray(edges, dtype=float)
    lo, hi = edges[:-1], edges[1:]
    right_share = np.clip((hi - origin) / np.maximum(hi - lo, EPS), 0.0, 1.0)
    right = np.sum(p * right_share, axis=-1)
    return 2.0 * right - np.sum(p, axis=-1)


def wasserstein1(p: np.ndarray, q: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Exact ∫|F_p − F_q| dx for piecewise-linear CDFs on shared edges (bp)."""
    w = widths(edges)
    d = cdf_at_edges(p) - cdf_at_edges(q)
    d0, d1 = d[..., :-1], d[..., 1:]
    same = d0 * d1 >= 0
    a0, a1 = np.abs(d0), np.abs(d1)
    with np.errstate(invalid="ignore", divide="ignore"):
        crossing = (d0**2 + d1**2) / (2 * np.maximum(a0 + a1, EPS))
    per_bin = np.where(same, (a0 + a1) / 2, crossing) * w
    return np.sum(per_bin, axis=-1)


def wasserstein2(p: np.ndarray, q: np.ndarray, edges: np.ndarray, grid: int = 2000) -> np.ndarray:
    """sqrt ∫₀¹ (Q_p(u) − Q_q(u))² du, evaluated on a midpoint grid of u (bp)."""
    u = (np.arange(grid) + 0.5) / grid
    edges = np.asarray(edges, dtype=float)
    p2 = np.asarray(p, dtype=float).reshape(-1, np.shape(p)[-1])
    q2 = np.asarray(q, dtype=float).reshape(-1, np.shape(q)[-1])
    out = np.empty(p2.shape[0])
    for i, (a, b) in enumerate(zip(p2, q2, strict=True)):
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
            out[i] = np.nan
            continue
        qa = _quantile_curve(a, edges, u)
        qb = _quantile_curve(b, edges, u)
        out[i] = np.sqrt(np.mean((qa - qb) ** 2))
    return out.reshape(np.shape(p)[:-1])


def _quantile_curve(p, edges, u):
    F = cdf_at_edges(p)
    keep = np.concatenate([[True], np.diff(F) > EPS])  # drop flat stretches for interpolation
    return np.interp(u, F[keep], edges[keep])


def jensen_shannon(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Jensen–Shannon divergence in bits (0 = identical, 1 = disjoint)."""
    m = (p + q) / 2

    def kl(a, b):
        with np.errstate(invalid="ignore", divide="ignore"):
            terms = np.where(
                a > 0, a * np.log2(np.where(a > 0, a, 1.0) / np.where(b > 0, b, 1.0)), 0.0
            )
        return terms.sum(axis=-1)

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


# ---------------------------------------------------------------- background & amplitude


def flank_mask(edges: np.ndarray, width_bp: float) -> np.ndarray:
    """Bins lying entirely within `width_bp` of either window end."""
    edges = np.asarray(edges, dtype=float)
    lo, hi = edges[0], edges[-1]
    return (edges[1:] <= lo + width_bp + EPS) | (edges[:-1] >= hi - width_bp - EPS)


def core_mask(edges: np.ndarray, half_width_bp: float) -> np.ndarray:
    """Bins lying entirely within ±half_width of the anchor."""
    edges = np.asarray(edges, dtype=float)
    return (edges[:-1] >= -half_width_bp - EPS) & (edges[1:] <= half_width_bp + EPS)


def background(
    values: np.ndarray, edges: np.ndarray, width_bp: float, estimator: str = "mean"
) -> np.ndarray:
    """Per-profile background from the outer flanks (NaN-aware).

    ``mean`` is unbiased for the background level of skewed count noise; ``median`` is robust to
    a neighbouring peak in the flank but sits below the mean of right-skewed noise, leaving a
    positive floor after subtraction.
    """
    mask = flank_mask(edges, width_bp)
    if not mask.any():
        raise ValueError("Background flank width is smaller than one bin")
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        flank = values[..., mask]
        if estimator == "median":
            return np.nanmedian(flank, axis=-1)
        if estimator == "mean":
            return np.nanmean(flank, axis=-1)
        raise ValueError("background estimator must be mean or median")


def amplitude_table(
    values: np.ndarray, edges: np.ndarray, center_half: float, flank_width: float
) -> dict:
    """Amplitude summaries for (..., bins) signal; all integrals are signal × bp."""
    w = widths(edges)
    finite = np.isfinite(values)
    x = np.where(finite, values, 0.0)
    covered_bp = np.sum(finite * w, axis=-1)
    auc = np.sum(x * w, axis=-1)
    centre = core_mask(edges, center_half)
    flank = flank_mask(edges, flank_width)
    centre_auc = np.sum((x * w)[..., centre], axis=-1)
    flank_auc = np.sum((x * w)[..., flank], axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        centre_mean = centre_auc / np.sum(w[centre])
        flank_mean = flank_auc / np.sum(w[flank])
        return dict(
            auc=auc,
            mean=np.where(covered_bp > 0, auc / np.maximum(covered_bp, EPS), np.nan),
            max=np.where(
                finite.any(axis=-1), np.max(np.where(finite, values, -np.inf), axis=-1), np.nan
            ),
            centre_auc=centre_auc,
            flank_auc=flank_auc,
            centre_flank_ratio=(centre_mean + EPS) / (flank_mean + EPS),
            missing_fraction=1 - finite.mean(axis=-1),
            zero_fraction=np.mean(finite & (values == 0), axis=-1),
            negative_fraction=np.mean(finite & (values < 0), axis=-1),
        )


def shape_table(p: np.ndarray, edges: np.ndarray, oriented: bool) -> dict:
    """Per-profile shape metrics; signed metrics are NaN for unoriented anchors."""
    table = dict(
        centroid=centroid(p, edges),
        sd_width=spread_sd(p, edges),
        iqr_width=iqr_width(p, edges),
        asymmetry=asymmetry(p, edges),
    )
    if not oriented:
        table["centroid"] = np.full_like(table["centroid"], np.nan)
        table["asymmetry"] = np.full_like(table["asymmetry"], np.nan)
    return table
