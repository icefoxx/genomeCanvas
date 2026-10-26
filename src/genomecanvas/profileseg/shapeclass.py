"""Condition-blind shape classes and anchor alignment (ideas from FunChIP, Parodi et al. 2017).

FunChIP clusters ChIP-seq peaks by shape with a shift-only k-mean alignment (k-mean alignment,
as in the fdakma package). Here the same idea serves inference, always on the
**pooled** profile of each region (all samples together, condition labels unused), so that neither
the classes nor the alignment can absorb or create a condition difference:

* shape classes give an optional descriptive table of which kinds of profiles change (using them to
  stratify the pooled permutation null was evaluated and rejected: no calibration gain, less power);
* the per-region shift that aligns a region to its class template corrects anchor misplacement
  (summit / motif uncertainty), which otherwise blurs aggregate profiles. The same shift is applied
  to every sample of the region.

Distances are L2 between mass-normalised, non-negative profiles on the core window, optionally
plus a first-derivative term (FunChIP's d1). The number of classes is chosen by the simplified
silhouette (distance to own vs nearest other template); weak structure gives one class.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import metrics as M
from .data import Profiles


@dataclass
class ShapeClasses:
    labels: np.ndarray        # region → class (0..k-1), for the regions passed in
    shift_bins: np.ndarray    # integer shift aligning each region to its template (pooled profile)
    templates: np.ndarray     # k × bins
    k: int
    silhouette: float
    silhouettes: dict         # k → simplified silhouette


def _shifted(P: np.ndarray, q: int) -> np.ndarray:
    """Shift rows of P by q bins (positive = move signal to higher positions), zero fill."""
    out = np.zeros_like(P)
    if q > 0:
        out[:, q:] = P[:, :-q]
    elif q < 0:
        out[:, :q] = P[:, -q:]
    else:
        out[:] = P
    return out


def _best_shift(P: np.ndarray, T: np.ndarray, Q: int, alpha: float = 0.0, w: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """For each row of P, the shift in [-Q, Q] minimising the distance to template T."""
    best_d = np.full(P.shape[0], np.inf)
    best_q = np.zeros(P.shape[0], int)
    dT = np.gradient(T)
    for q in range(-Q, Q + 1):
        S = _shifted(P, q)
        d = ((S - T) ** 2).sum(1)
        if alpha > 0:
            d = (1 - alpha) * d + alpha * w * ((np.gradient(S, axis=1) - dT) ** 2).sum(1)
        better = d < best_d
        best_d[better], best_q[better] = d[better], q
    return best_q, best_d


def _align_rows(P: np.ndarray, shifts: np.ndarray) -> np.ndarray:
    out = np.empty_like(P)
    for q in np.unique(shifts):
        m = shifts == q
        out[m] = _shifted(P[m], int(q))
    return out


def kmean_alignment(P: np.ndarray, k: int, max_shift: int, seed: int = 0, iterations: int = 15,
                    alpha: float = 0.0, w: float = 1.0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Shift-only k-mean alignment of the rows of P (regions × bins, mass-normalised)."""
    rng = np.random.default_rng(seed)
    n = P.shape[0]
    # k-means++ initialisation on unaligned profiles
    centres = [P[rng.integers(n)]]
    for _ in range(1, k):
        d = np.min([((P - c) ** 2).sum(1) for c in centres], axis=0)
        centres.append(P[rng.choice(n, p=d / d.sum())] if d.sum() > 0 else P[rng.integers(n)])
    T = np.array(centres)
    labels = np.zeros(n, int)
    shifts = np.zeros(n, int)
    for _ in range(iterations):
        D = np.empty((k, n))
        Qs = np.empty((k, n), int)
        for j in range(k):
            Qs[j], D[j] = _best_shift(P, T[j], max_shift, alpha, w)
        new = D.argmin(0)
        shifts = Qs[new, np.arange(n)]
        changed = (new != labels).mean()
        labels = new
        for j in range(k):
            m = labels == j
            if not m.any():
                T[j] = P[rng.integers(n)]
                continue
            T[j] = _align_rows(P[m], shifts[m]).mean(0)
            # average warping = identity: re-centre the template on the members' mean shift
            mean_q = int(round(shifts[m].mean()))
            if mean_q:
                T[j] = _shifted(T[j][None], -mean_q)[0]
                shifts[m] -= mean_q
        if changed < 1e-3:
            break
    D = np.array([_best_shift(P, T[j], max_shift, alpha, w)[1] for j in range(k)])
    return labels, shifts, T, D


def simplified_silhouette(D: np.ndarray, labels: np.ndarray) -> float:
    if D.shape[0] < 2:
        return 0.0
    n = D.shape[1]
    a = D[labels, np.arange(n)]
    Dm = D.copy()
    Dm[labels, np.arange(n)] = np.inf
    b = Dm.min(0)
    s = (b - a) / np.maximum(np.maximum(a, b), 1e-12)
    return float(np.mean(s))


def shape_classes(P: np.ndarray, max_shift: int, k_max: int = 6, min_size: int = 30, seed: int = 0,
                  min_silhouette: float = 0.15, alpha: float = 0.0) -> ShapeClasses:
    """Choose k in 2..k_max by simplified silhouette; one class when structure is weak."""
    P = np.asarray(P, float)
    w = 1.0
    if alpha > 0:  # FunChIP: derivative weight = median ratio of the two distances
        d0 = ((P - P.mean(0)) ** 2).sum(1)
        d1 = ((np.gradient(P, axis=1) - np.gradient(P.mean(0))) ** 2).sum(1)
        w = float(np.median(d0) / max(np.median(d1), 1e-12))
    lab1, sh1, T1, D1 = kmean_alignment(P, 1, max_shift, seed, alpha=alpha, w=w)
    best = ShapeClasses(lab1, sh1, T1, 1, 0.0, {1: 0.0})
    sil = {1: 0.0}
    for k in range(2, k_max + 1):
        if P.shape[0] < k * min_size:
            break
        lab, sh, T, D = kmean_alignment(P, k, max_shift, seed + k, alpha=alpha, w=w)
        sizes = np.bincount(lab, minlength=k)
        if sizes.min() < min_size:
            sil[k] = float("nan")
            continue
        s = simplified_silhouette(D, lab)
        sil[k] = s
        if s > best.silhouette and s >= min_silhouette:
            best = ShapeClasses(lab, sh, T, k, s, sil)
    best.silhouettes = sil
    if best.k == 1:  # alignment shifts from the single-template fit are still useful
        best.shift_bins = sh1
    return best


def pooled_core(values: np.ndarray, edges: np.ndarray, background_width: int, core_flank: int,
                estimator: str = "mean") -> tuple[np.ndarray, np.ndarray]:
    """Mass-normalised pooled (all samples) background-corrected core profiles, and a usable mask."""
    bg = M.background(values, edges, background_width, estimator)
    core = M.core_mask(edges, core_flank)
    y = np.clip(np.nan_to_num((values - bg[..., None])[..., core]).mean(0), 0, None)
    mass = y.sum(-1)
    ok = mass > 0
    P = np.zeros_like(y)
    P[ok] = y[ok] / mass[ok, None]
    return P, ok


def align_profiles(profiles: Profiles, background_width: int, core_flank: int, max_shift_bp: float = 150.0,
                   seed: int = 0, estimator: str = "mean") -> tuple[Profiles, np.ndarray, ShapeClasses | None]:
    """Re-centre each region on its class template (shift from the pooled profile, all samples alike).

    Returns the aligned profiles, the shift applied per region in bp (0 for regions without
    positive pooled signal), and the shape classes used.
    """
    edges = profiles.edges
    width = float(np.median(np.diff(edges)))
    Q = max(1, int(round(max_shift_bp / width)))
    P, ok = pooled_core(profiles.values, edges, background_width, core_flank, estimator)
    shift_bp = np.zeros(len(profiles.regions))
    if ok.sum() < 30:
        return profiles, shift_bp, None
    sc = shape_classes(P[ok], Q, seed=seed)
    q = np.zeros(len(profiles.regions), int)
    q[ok] = sc.shift_bins
    v = profiles.values.copy()
    for r in np.flatnonzero(q):
        s = q[r]
        row = v[:, r]
        out = np.empty_like(row)
        if s > 0:  # signal moves right: pad the left edge with the first value
            out[:, s:] = row[:, :-s]
            out[:, :s] = row[:, :1]
        else:
            out[:, :s] = row[:, -s:]
            out[:, s:] = row[:, -1:]
        v[:, r] = out
    shift_bp = q * width
    aligned = Profiles(v, edges, profiles.samples, profiles.regions, profiles.oriented,
                       provenance=dict(profiles.provenance, aligned_max_shift_bp=max_shift_bp))
    return aligned, shift_bp, sc
