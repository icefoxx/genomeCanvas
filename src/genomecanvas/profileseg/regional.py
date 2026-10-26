"""ProfileSeg v2: per-region tests with empirical-Bayes moderated statistics.

Each region is tested with biological replicates as observations, using two-group linear models
whose residual variances are shrunk towards a common prior estimated across regions (Smyth 2004,
*Stat Appl Genet Mol Biol* 3:3; the method behind limma). Borrowing information across regions is
what makes per-region inference feasible with 2–3 replicates. Multiplicity across regions is
controlled with Benjamini–Hochberg FDR, separately for each metric.

This answers "which regions change?". It complements, and does not replace, the global
replicate-level summary in ``analysis.py``.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy import special, stats

from . import metrics as M
from .analysis import Analysis, _welch

REGION_METRICS = {
    "log2_auc": "absolute amplitude: log2(window AUC); depends on upstream normalisation",
    "log2_signal": "log2 background-corrected core signal (absolute; depends on upstream normalisation)",
    "log2_signal_relative": "log2 background-corrected core signal relative to the bulk of regions (site-specific; assumes most regions unchanged)",
    "log2_enrichment": "log2(centre signal / local background) per region; insensitive to global rescaling",
    "shape_kernel": "omnibus shape test (kernel distance between condition-mean profiles; pooled permutation null); no direction",
    "shift": "directional shift test (projection on the template derivative, pooled permutation null); estimate = shift in bp from a template fit",
    "centroid": "centre of mass of the background-corrected core profile (bp; oriented anchors); low power on broad profiles",
    "mad_width": "spread (mean absolute deviation, bp) of the background-corrected core profile; low power, exploratory",
}

# Tuned on separate seeds (benchmarks/profileseg_region_tuning.py); see docs/PROFILESEG_COMPARISON.md §3.
KERNEL_SUPPORT_FRACTION = None  # full core window: restricting to the support lost power on real shapes
KERNEL_SIGMA_BP = 60.0
KERNEL_S0_QUANTILE = 0.5
TEMPLATE_SMOOTH = 5  # template of the shift score test (calibration validated with this value)
# Depth diagnostic: per-region shape tests had power with ~400 fragments per region and sample (yeast
# MNase) and none with ~10 (TF CUT&RUN); below this median the report warns.
LOW_DEPTH_FRAGMENTS = 50
# Effect-size estimator for `shift`: stronger smoothing and condition-blind template refinement
# (benchmarks/profileseg_shift_estimate_tuning.py: worst-cell estimate/true ratio 0.58 → 0.73).
SHIFT_FIT_SMOOTH = 7
SHIFT_FIT_REFINE = 3


# ---------------------------------------------------------------- empirical Bayes (Smyth 2004)


def trigamma_inverse(x: np.ndarray) -> np.ndarray:
    """Solve trigamma(y) = x (Newton iteration as in limma::trigammaInverse)."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    y = np.empty_like(x)
    big, small = x > 1e7, x < 1e-6
    y[big] = 1 / np.sqrt(x[big])
    y[small] = 1 / x[small]
    mid = ~(big | small)
    if mid.any():
        xm = x[mid]
        ym = 0.5 + 1 / xm
        for _ in range(50):
            tri = special.polygamma(1, ym)
            dif = tri * (1 - tri / xm) / special.polygamma(2, ym)
            ym = ym + dif
            if np.max(-dif / ym) < 1e-8:
                break
        y[mid] = ym
    return y


def _smooth_trend(x: np.ndarray, z: np.ndarray, bins: int = 10) -> np.ndarray:
    """Smooth trend of z on x: medians in quantile bins of x, linearly interpolated (flat at the ends)."""
    if x.size < 4 * bins:
        bins = max(1, x.size // 4)
    q = np.unique(np.quantile(x, np.linspace(0, 1, bins + 1)))
    if q.size < 3:
        return np.full_like(z, np.mean(z))
    which = np.clip(np.searchsorted(q, x, side="right") - 1, 0, q.size - 2)
    cx = np.array([np.median(x[which == b]) for b in range(q.size - 1)])
    cz = np.array([np.mean(z[which == b]) for b in range(q.size - 1)])
    return np.interp(x, cx, cz)


def fit_f_dist(s2: np.ndarray, df: float, covariate: np.ndarray | None = None) -> tuple[float, np.ndarray | float]:
    """Prior degrees of freedom d0 and scale s0² for sample variances s2 with df degrees of freedom.

    Follows limma::fitFDist. With ``covariate`` the prior scale follows a smooth trend in the
    covariate (as limma's ``trend=TRUE``, here with a binned smoother instead of a spline), which is
    needed when variance depends on signal strength.
    """
    s2 = np.asarray(s2, dtype=float)
    ok = np.isfinite(s2) & (s2 > 0)
    if ok.sum() < 3:
        return 0.0, float(np.nanmedian(s2)) if np.isfinite(s2).any() else 1.0
    floor = 1e-5 * np.median(s2[ok])
    z = np.log(np.maximum(s2[ok], floor))
    e = z - special.digamma(df / 2) + np.log(df / 2)
    if covariate is None:
        trend = np.full_like(e, np.mean(e))
    else:
        trend = _smooth_trend(np.asarray(covariate, dtype=float)[ok], e)
    resid = e - trend
    evar = float(np.sum(resid**2) / (resid.size - 1)) - float(special.polygamma(1, df / 2))
    if evar > 0:
        d0 = float(2 * trigamma_inverse(evar)[0])
        scale = np.exp(trend + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0, scale = np.inf, np.exp(trend)
    if covariate is None:
        return d0, float(scale[0])
    # The scale is a function of the covariate; interpolating also covers zero / missing variances.
    cov = np.asarray(covariate, dtype=float)
    order = np.argsort(cov[ok])
    full = np.interp(cov, cov[ok][order], scale[order])
    return d0, full


def moderated_t(values: np.ndarray, ref: np.ndarray, alt: np.ndarray, covariate: np.ndarray | None = None) -> dict:
    """Two-group moderated t for every column of values (samples × regions)."""
    x, y = values[ref], values[alt]
    n1, n2 = len(ref), len(alt)
    diff = y.mean(axis=0) - x.mean(axis=0)
    d = n1 + n2 - 2
    if d < 1:
        nan = np.full(values.shape[1], np.nan)
        return dict(estimate=diff, t=nan, df=nan, p=nan, d0=np.nan, s02=np.nan)
    s2 = (((x - x.mean(axis=0)) ** 2).sum(axis=0) + ((y - y.mean(axis=0)) ** 2).sum(axis=0)) / d
    d0, s02 = fit_f_dist(s2, d, covariate)
    if np.isinf(d0):
        post, df_total = np.broadcast_to(s02, s2.shape).astype(float), np.full_like(s2, np.inf)
    else:
        post = (d0 * s02 + d * s2) / (d0 + d)
        df_total = np.full_like(s2, d0 + d)
    se = np.sqrt(post * (1 / n1 + 1 / n2))
    t = diff / se
    p = 2 * (stats.norm.sf(np.abs(t)) if np.isinf(d0) else stats.t.sf(np.abs(t), df_total))
    return dict(estimate=diff, t=t, df=df_total, p=p, d0=d0, s02=s02 if np.ndim(s02) == 0 else float(np.nanmedian(s02)))


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg adjusted p-values (NaN preserved)."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if not ok.size:
        return out
    order = ok[np.argsort(p[ok])]
    m = len(order)
    adj = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out[order] = np.minimum(adj, 1.0)
    return out


# ---------------------------------------------------------------- per-region metrics


@dataclass
class RegionalResult:
    regions: np.ndarray  # indices into profiles.regions
    tests: dict  # metric → dict(estimate, t, df, p, fdr, call, region_index)
    robust_global: list = field(default_factory=list)
    shape_regions: int = 0
    notes: list = field(default_factory=list)
    strength_trend: float = float("nan")  # Spearman(relative change, signal strength)
    influence: dict = field(default_factory=dict)  # replicate robustness diagnostics (see replicate_influence)
    classes: list = field(default_factory=list)  # condition-blind shape classes: per-class summaries
    class_templates: np.ndarray | None = None
    class_labels: np.ndarray | None = None  # per shape-tested region (order of the shape tests)
    depth: dict = field(default_factory=dict)  # fragments per region and sample (exact for BAM, estimated for bigWig)


def amplitude_matrices(a: Analysis) -> tuple[np.ndarray, dict]:
    """Per sample × region amplitude metrics for regions complete in every sample."""
    c = a.config
    comp = np.flatnonzero(a.complete)
    auc = a.amplitude["auc"][:, comp]
    pc = 1e-3 * max(float(np.median(auc[auc > 0])) if np.any(auc > 0) else 1.0, M.EPS)
    centre = a.amplitude["centre_auc"][:, comp]
    bgw = a.background[:, comp] * (2 * c.center_half)
    # Background-corrected signal: a uniform change in IP efficiency scales it equally in every
    # region, whereas window AUC (which contains background) changes less in weak regions.
    mass = a.core_mass[:, comp]
    positive = mass[mass > 0]
    pc_signal = 0.1 * float(np.median(positive)) if positive.size else 1.0
    return comp, dict(
        log2_auc=np.log2(np.maximum(auc, 0) + pc),
        log2_signal=np.log2(np.maximum(mass, 0) + pc_signal),
        log2_enrichment=np.log2(np.maximum(centre, 0) + pc) - np.log2(np.maximum(bgw, 0) + pc),
    )


def support_mask(y: np.ndarray, fraction: float = 0.2) -> np.ndarray:
    """Condition-blind support per region: bins where the pooled profile exceeds ``fraction`` of its max.

    ``y`` is samples × regions × bins (background-corrected). Using all samples together keeps the
    region definition independent of the contrast.
    """
    pooled = np.nanmean(y, axis=0)
    return pooled > fraction * np.nanmax(pooled, axis=1, keepdims=True)


def support_moments(y: np.ndarray, edges: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centroid and mean absolute deviation (bp) of the positive signal inside the support.

    Moments over the full window are dominated by noise far from the anchor (its leverage grows
    with distance); restricting to the support removes most of it. The restriction attenuates
    effect sizes, so these are test statistics; effect sizes come from :func:`template_fit`.
    """
    p = M.normalise(np.clip(np.nan_to_num(y), 0, None) * mask * M.widths(edges))
    return M.centroid(p, edges), M.mean_abs_deviation(p, edges)


def _interp_rows(table: np.ndarray, pos: np.ndarray, x0: float, h: float) -> np.ndarray:
    """Linear interpolation of each row of ``table`` (regions × bins on a uniform grid) at ``pos``
    (… × regions × bins, in bp); zero outside the grid."""
    f = (pos - x0) / h
    i = np.floor(f).astype(int)
    w = f - i
    nb = table.shape[-1]
    lo_ok, hi_ok = (i >= 0) & (i < nb), (i + 1 >= 0) & (i + 1 < nb)
    t = np.broadcast_to(table, pos.shape)
    lo = np.where(lo_ok, np.take_along_axis(t, np.clip(i, 0, nb - 1), -1), 0.0)
    hi = np.where(hi_ok, np.take_along_axis(t, np.clip(i + 1, 0, nb - 1), -1), 0.0)
    return lo * (1 - w) + hi * w


def template_fit(y: np.ndarray, cov: np.ndarray, edges: np.ndarray, smooth: int = 5, iterations: int = 4,
                 floor: float = 0.05, chunk: int = 2000, refine: int = 0) -> dict:
    """Per-sample shift and log2 scale relative to each region's pooled template (see
    :func:`_template_fit_block`); regions are processed in blocks to bound memory."""
    R = y.shape[1]
    if R <= chunk:
        return _template_fit_block(y, cov, edges, smooth, iterations, floor, refine)
    parts = [_template_fit_block(y[:, i:i + chunk], cov[:, i:i + chunk], edges, smooth, iterations, floor, refine)
             for i in range(0, R, chunk)]
    return {k: np.concatenate([q[k] for q in parts], axis=1) for k in parts[0]}


def _template_fit_block(y: np.ndarray, cov: np.ndarray, edges: np.ndarray, smooth: int = 5, iterations: int = 4,
                        floor: float = 0.05, refine: int = 0) -> dict:
    """Per-sample shift (bp) and log2 scale of each region relative to its pooled template.

    Model: y_s(x) ≈ k_s · T((x − μ − δ_s) / g_s + μ), with T the condition-blind pooled profile
    (lightly smoothed) and μ its centroid. Weighted Gauss–Newton (weights 1/coverage, Poisson-like)
    on bins where T exceeds ``floor`` of its maximum. Differences of δ between conditions estimate
    the shift with little attenuation, unlike support-restricted moments.

    The pooled template of a region whose conditions differ is a blurred mixture, which attenuates
    the estimated shifts. With ``refine`` > 0 the template is rebuilt from the samples re-aligned by
    their own estimated shifts (no condition labels used) and the fit repeated.
    """
    c = M.centers(edges)
    h = float(np.median(np.diff(c)))
    y = np.nan_to_num(y)
    S, R, _ = y.shape
    kern = None
    if smooth > 1:
        half = (smooth + 1) // 2
        kern = np.r_[np.arange(1, half + 1), np.arange(half - 1, 0, -1)].astype(float)
        kern /= kern.sum()

    def make_template(profiles):
        T = np.clip(profiles.mean(axis=0), 0, None)
        return T if kern is None else np.apply_along_axis(lambda r: np.convolve(r, kern, mode="same"), 1, T)

    T0 = make_template(y)
    delta, logg = np.zeros((S, R)), np.zeros((S, R))
    for stage in range(refine + 1):
        if stage:  # re-align each sample by its own shift (relative to the mean shift) and rebuild T
            back = delta - delta.mean(0, keepdims=True)
            aligned = _interp_rows(y.reshape(S * R, -1), (c + back.reshape(-1, 1)), c[0], h).reshape(y.shape)
            T0 = make_template(aligned)
            delta, logg = np.zeros((S, R)), np.zeros((S, R))
        mu = (T0 * c).sum(1) / np.maximum(T0.sum(1), 1e-12)
        use = (T0 > floor * T0.max(1, keepdims=True))[None]
        wt = np.where(use, 1.0 / np.maximum(np.nan_to_num(cov), 1.0), 0.0)
        mu_ = mu[None, :, None]
        for _ in range(iterations):
            g = 2.0 ** logg[..., None]
            xs = (c - mu_ - delta[..., None]) / g + mu_
            T = _interp_rows(T0, xs, c[0], h)
            dT = np.gradient(T, h, axis=-1)
            X = np.stack([T, -dT / g, -(xs - mu_) * dT * np.log(2)], axis=-1)  # S × R × B × 3
            XtW = X * wt[..., None]
            A = np.einsum("srbi,srbj->srij", XtW, X) + 1e-9 * np.eye(3)
            b = np.einsum("srbi,srb->sri", XtW, y)
            beta = np.linalg.solve(A, b[..., None])[..., 0]
            k = np.where(np.abs(beta[..., 0]) > 1e-9, beta[..., 0], 1e-9)
            delta = np.clip(delta + beta[..., 1] / k, -400, 400)
            logg = np.clip(logg + beta[..., 2] / k, -1.5, 1.5)
    return dict(shift=delta, log2_scale=logg)


def label_assignments(n1: int, n: int, blocks: np.ndarray | None = None) -> tuple[list[tuple], int]:
    """Distinct replicate-label assignments (sets of indices called group A) and the index of the
    observed one (the first n1 samples). An assignment and its mirror count once when n1 = n − n1.

    With ``blocks`` (processing batches, one label per sample, in the same order) only assignments
    that keep the observed number of group-A samples in every block are allowed: labels are
    permuted within batches, so batch differences stay out of the null. If only the observed
    assignment remains (e.g. condition and batch coincide), the test cannot be performed.
    """
    splits = [s for s in itertools.combinations(range(n), n1) if n1 != n - n1 or 0 in s]
    if blocks is not None:
        blocks = np.asarray(blocks)
        obs = {b: int(np.sum(blocks[:n1] == b)) for b in np.unique(blocks)}
        def keeps(sel):
            return all(int(np.sum(blocks[list(sel)] == b)) == k for b, k in obs.items())
        splits = [s for s in splits
                  if keeps(s) or (n1 == n - n1 and keeps(tuple(np.setdiff1d(np.arange(n), s))))]
    return splits, splits.index(tuple(range(n1)))


def pooled_permutation_p(F: np.ndarray, obs_i: int, strata: np.ndarray | None = None, min_size: int = 30) -> np.ndarray:
    """p-values of the observed row of F (assignments × regions) against the non-observed rows,
    pooled over regions. ``strata`` (experimental, not used by ``per_region_tests``) pools within
    strata of at least ``min_size`` regions; it was evaluated as shape-class stratification and
    rejected (docs/PROFILESEG_COMPARISON.md §7)."""
    obs = F[obs_i]
    rest = np.delete(F, obs_i, axis=0)
    p = np.empty(obs.shape)
    groups = np.zeros(obs.shape, int) if strata is None else np.asarray(strata)
    small = np.isin(groups, [g for g in np.unique(groups) if np.sum(groups == g) < min_size])
    for g in np.unique(groups[~small]):
        m = (groups == g) & ~small
        null = np.sort(rest[:, m].ravel())
        p[m] = (1 + null.size - np.searchsorted(null, obs[m], side="left")) / (1 + null.size)
    if small.any():
        null = np.sort(rest.ravel())
        p[small] = (1 + null.size - np.searchsorted(null, obs[small], side="left")) / (1 + null.size)
    return p


def kernel_shape_test(y: np.ndarray, mask: np.ndarray, ref: np.ndarray, alt: np.ndarray, sigma_bp: float = 60.0,
                      s0_quantile: float = 0.5, centers: np.ndarray | None = None,
                      strata: np.ndarray | None = None, blocks: np.ndarray | None = None) -> np.ndarray:
    """Omnibus per-region shape test: a kernel (MMD-type) distance between condition means of the
    mass-normalised profiles, scaled by within-condition spread, with a null pooled over regions
    and all non-observed replicate-label assignments.

    One region alone has too few label assignments (3 vs 3: 10; 2 vs 2: 3); pooling the
    assignments over regions gives a usable null, as in SAM-type permutation schemes. The Gaussian
    kernel is local, so noise far from the signal has little influence.

    The stabilising constant s0 must not depend on the label assignment: assignments that mix the
    conditions carry part of a true effect in their within-group spread, so an s0 taken per
    assignment is larger for the null assignments than for the observed one and makes the test
    anti-conservative once many regions change. s0 is therefore a quantile of the total spread
    (identical under every assignment), scaled to a per-sample variance.
    """
    # Normalisation: the signed (unclipped) background-corrected profile divided by its mass inside the
    # condition-blind support of the pooled profile. Clipping noise at zero added a positive floor whose
    # relative size grows as the signal weakens, so amplitude changes looked like broadening (4 vs 4,
    # 10 bp bins: 1–3 % of regions called under amplitude-only changes). Dividing by the whole-window
    # signed mass removed that bias but the near-zero, noisy denominator of broad profiles destroyed
    # power; the support mass is a stable denominator that still cancels pure amplitude changes.
    # Regions whose support mass is not positive in every sample are not tested (NaN).
    y = np.nan_to_num(y) * mask
    support = support_mask(y, 0.2)
    mass = (y * support[None]).sum(-1)
    testable = np.all(mass > 0, axis=0)
    p = np.where(testable[None, :, None], y / np.where(mass > 0, mass, 1.0)[..., None], 0.0)
    c = np.arange(p.shape[-1], dtype=float) if centers is None else centers
    K = np.exp(-0.5 * ((c[:, None] - c[None, :]) / sigma_bp) ** 2)
    # Kernel feature map: K = L Lᵀ (eigendecomposition; the Gaussian kernel is positive semi-definite),
    # so xᵀ K x = ‖Lᵀ x‖². Distances become Euclidean in φ = p L: identical statistics, and each label
    # assignment costs O(regions × bins) instead of O(regions × bins²).
    lam, V = np.linalg.eigh(K)
    keep = lam > 1e-10 * lam.max()  # numerically zero directions contribute nothing; a smooth kernel has low rank
    L = V[:, keep] * np.sqrt(lam[keep])
    n1 = len(ref)
    n = n1 + len(alt)
    phi = p[np.r_[ref, alt]][:, testable] @ L  # untestable regions take no part in s0 or the null
    total = ((phi - phi.mean(0)) ** 2).sum((0, 2))
    s0 = np.quantile(total / (n - 1), s0_quantile)
    splits, obs_i = label_assignments(n1, n, None if blocks is None else np.asarray(blocks)[np.r_[ref, alt]])
    if len(splits) < 2 or not testable.any():
        return np.full(p.shape[1], np.nan)
    F = []
    for s in splits:
        a = np.array(s)
        b = np.setdiff1d(np.arange(n), a)
        ma, mb = phi[a].mean(0), phi[b].mean(0)
        D = ((ma - mb) ** 2).sum(-1)
        W = (((phi[a] - ma) ** 2).sum((0, 2)) + ((phi[b] - mb) ** 2).sum((0, 2))) / max(n - 2, 1)
        F.append(D / (W + s0))
    F = np.array(F)
    pval = np.full(p.shape[1], np.nan)
    if testable.any():
        pval[testable] = pooled_permutation_p(F, obs_i, None if strata is None else np.asarray(strata)[testable])
    return pval


def shift_score_test(y: np.ndarray, ref: np.ndarray, alt: np.ndarray, centers: np.ndarray, smooth: int = 5,
                     s0_quantile: float = 0.5, strata: np.ndarray | None = None,
                     blocks: np.ndarray | None = None) -> dict:
    """Directional per-region shift test: each sample's mass-normalised profile is projected on the
    derivative of the pooled (condition-blind, lightly smoothed) template, the locally most powerful
    direction for a small translation. The group difference of the projections is scaled by their
    within-condition variance plus a label-invariant s0, with the same pooled permutation null as
    :func:`kernel_shape_test`.

    For 50 bp shifts this is 2–4 times as powerful as the omnibus kernel test (it spends one degree
    of freedom instead of many); it does not respond to symmetric width changes. Returns p-values and
    the sign of the shift (alternative − reference). The projection is attenuated for shifts that are
    large relative to the peak, so effect sizes come from :func:`template_fit`.
    """
    yc = np.clip(np.nan_to_num(y), 0, None)
    p = yc / np.maximum(yc.sum(-1, keepdims=True), 1e-12)
    p = p[np.r_[ref, alt]]
    n1, n = len(ref), len(ref) + len(alt)
    Tm = p.mean(0)
    if smooth > 1:
        half = (smooth + 1) // 2
        kern = np.r_[np.arange(1, half + 1), np.arange(half - 1, 0, -1)].astype(float)
        kern /= kern.sum()
        Tm = np.apply_along_axis(lambda r: np.convolve(r, kern, mode="same"), 1, Tm)
    g = np.gradient(Tm, float(np.median(np.diff(centers))), axis=-1)
    z = -((p - Tm) * g).sum(-1) / np.maximum((g * g).sum(-1), 1e-12)  # ≈ per-sample shift (bp), attenuated
    s0 = np.quantile(((z - z.mean(0)) ** 2).sum(0) / (n - 1), s0_quantile)
    splits, obs_i = label_assignments(n1, n, None if blocks is None else np.asarray(blocks)[np.r_[ref, alt]])
    if len(splits) < 2:
        return dict(p=np.full(z.shape[1], np.nan), sign=np.sign(z[n1:].mean(0) - z[:n1].mean(0)))
    F = []
    for sp in splits:
        a = np.array(sp)
        b = np.setdiff1d(np.arange(n), a)
        d = z[b].mean(0) - z[a].mean(0)
        W = (((z[a] - z[a].mean(0)) ** 2).sum(0) + ((z[b] - z[b].mean(0)) ** 2).sum(0)) / max(n - 2, 1)
        F.append(d * d / (W + s0))
    F = np.array(F)
    pval = pooled_permutation_p(F, obs_i, strata)
    return dict(p=pval, sign=np.sign(z[n1:].mean(0) - z[:n1].mean(0)))


def replicate_influence(y: np.ndarray, cov: np.ndarray, edges: np.ndarray, ref: np.ndarray, alt: np.ndarray,
                        calls: dict, fdr: float, oriented: bool, names: list[str] | None = None) -> dict:
    """Sample-level robustness of the per-region shape calls.

    Per-region tests assume replicates are exchangeable. A single replicate with a sample-wide
    artefact (different fragment length, peak-calling offset, poor IP) touches every region at
    once. Two checks:

    * sample-wide shape offsets (label-free): each sample's median per-region shift and log2 scale
      relative to the pooled template; a sample far from the others is flagged.
    * leave-one-replicate-out: the share of shape calls that remain when each replicate is dropped.
      Calls that disappear when one particular replicate is removed (but survive the removal of
      any other) depend on that replicate.
    """
    names = names or [f"s{i}" for i in range(len(ref) + len(alt))]
    order = np.r_[ref, alt]
    fit = template_fit(y[order], cov[order], edges, smooth=TEMPLATE_SMOOTH)
    med_shift = np.median(fit["shift"], axis=1)
    med_scale = np.median(fit["log2_scale"], axis=1)
    rows, flags = [], []
    c = M.centers(edges)
    full = {m: set(np.flatnonzero(np.isin(v, ("Higher", "Lower", "Changed")))) for m, v in calls.items()}
    for j, i in enumerate(order):
        keep = np.delete(np.arange(len(order)), j)
        n_ref_k = int(np.sum(keep < len(ref)))
        n_alt_k = len(keep) - n_ref_k
        row = dict(sample=names[i], group="reference" if j < len(ref) else "alternative",
                   median_shift_bp=float(med_shift[j]), median_log2_scale=float(med_scale[j]))
        yk = y[order][keep]
        rk, ak = np.arange(n_ref_k), np.arange(n_ref_k, len(keep))
        if n_ref_k >= 1 and n_alt_k >= 1 and _n_assignments(n_ref_k, n_alt_k) >= 3 and len(keep) >= 3:
            loo = {}
            if "shape_kernel" in full:
                loo["shape_kernel"] = kernel_shape_test(yk, np.ones(yk.shape[1:], bool), rk, ak, KERNEL_SIGMA_BP,
                                                        KERNEL_S0_QUANTILE, c)
            if "shift" in full and oriented:
                loo["shift"] = shift_score_test(yk, rk, ak, c, smooth=TEMPLATE_SMOOTH)["p"]
            for m, pv in loo.items():
                kept = set(np.flatnonzero(bh(pv) < fdr))
                row[f"{m}_retained"] = len(kept & full[m]) / len(full[m]) if full[m] else float("nan")
        rows.append(row)
    # Sample-wide offsets are judged against the replicate's own group: a treatment moves the whole
    # group, an artefact moves one replicate away from its group-mates. (Comparing with all samples
    # flagged genuine condition effects.) With two replicates per group the culprit cannot be
    # identified, so the pair is reported.
    groups = [list(range(len(ref))), list(range(len(ref), len(order)))]
    for key, unit, floor in (("median_shift_bp", "bp", 10.0), ("median_log2_scale", "log2", 0.1)):
        v = np.array([r[key] for r in rows])
        dev = np.full(len(v), np.nan)
        for g in groups:
            for i in g:
                if len(g) >= 3:  # median of the whole group, so one deviant replicate does not shift it
                    dev[i] = v[i] - np.median(v[g])
                elif len(g) == 2:
                    dev[i] = v[i] - v[[j for j in g if j != i][0]]
        for g in groups:
            for i in g:
                if not np.isfinite(dev[i]):
                    continue
                others = np.abs(np.delete(dev, [j for j in g]))
                others = others[np.isfinite(others)]
                ref_spread = np.median(others) if others.size else 0.0
                if abs(dev[i]) > max(floor, 3 * ref_spread):
                    who = rows[i]["sample"] if len(g) > 2 else " / ".join(rows[j]["sample"] for j in g)
                    msg = (f"sample-wide {key.replace('median_', '')} offset within a group: {who} "
                           f"differs from its group-mates by {dev[i]:+.3g} {unit}: check fragment length, alignment "
                           "and IP quality (sample-wide artefacts are absorbed by the pooled null, but effect sizes "
                           "and the aggregate analysis can be affected)")
                    if msg not in flags:
                        flags.append(msg)
    # Leave-one-out flags need ≥ 3 replicates per group: with two, dropping one leaves 1 vs n and
    # far fewer label assignments, so retention mostly measures lost power (on real yeast MNase data
    # this flagged an unremarkable replicate). The retention table is still reported.
    for m in ("shape_kernel", "shift"):
        ret = [r.get(f"{m}_retained") for r in rows if np.isfinite(r.get(f"{m}_retained", np.nan))]
        if min(len(ref), len(alt)) >= 3 and len(ret) >= 3 and len(full.get(m, ())) >= 10:
            worst = int(np.argmin(ret))
            rest = np.delete(np.array(ret), worst)
            if ret[worst] < 0.5 and np.median(rest) >= 0.8:
                who = [r for r in rows if np.isfinite(r.get(f"{m}_retained", np.nan))][worst]["sample"]
                flags.append(f"{m}: only {ret[worst]:.0%} of calls remain without replicate {who} "
                             f"(others: median {np.median(rest):.0%}); the calls depend on that replicate")
    return dict(rows=rows, flags=flags)


def half_sample_mode(x: np.ndarray) -> float:
    """Robust mode (Bickel & Frühwirth 2006): repeatedly keep the densest half of the sorted data."""
    x = np.sort(np.asarray(x, dtype=float)[np.isfinite(x)])
    while x.size > 3:
        h = (x.size + 1) // 2
        k = int(np.argmin(x[h - 1:] - x[: x.size - h + 1]))
        x = x[k:k + h]
    if x.size == 3:
        return float(x[:2].mean() if x[1] - x[0] < x[2] - x[1] else x[1:].mean() if x[1] - x[0] > x[2] - x[1] else x[1])
    return float(x.mean()) if x.size else 0.0


def relative_offset(log_auc: np.ndarray, ref, alt, iterations: int = 5, exclude_q: float = 0.2):
    """Per-sample offsets estimated on regions that do not appear to change (TMM-like idea).

    A per-sample median is pulled by site-specific changes (if 20 % of regions gain, the median
    rises and unchanged regions look 'Lower'). Start from the half-sample mode of each sample's log
    ratio to the across-sample mean (the mode follows the unchanged majority), then iterate: take
    the median over regions with q ≥ ``exclude_q``, retest, and stop when the set is stable.
    """
    n = log_auc.shape[1]
    ratio = log_auc - log_auc.mean(axis=0, keepdims=True)
    offset = np.array([[half_sample_mode(r)] for r in ratio])
    keep = np.ones(n, bool)
    for _ in range(iterations):
        res = moderated_t(log_auc - offset, ref, alt)
        candidate = ~(bh(res["p"]) < exclude_q)
        if candidate.sum() < 0.5 * n or np.array_equal(candidate, keep):
            break
        keep = candidate
        offset = np.median(ratio[:, keep], axis=1, keepdims=True)
    return offset, keep


def per_region_tests(a: Analysis, fdr: float = 0.05, shape_quantile: float = 0.0) -> RegionalResult:
    """Moderated per-region tests for amplitude and shape metrics."""
    c = a.config
    n_ref = len(a.ref_idx)
    ref = np.arange(n_ref)
    alt = np.arange(n_ref, n_ref + len(a.alt_idx))
    notes = []
    if len(ref) + len(alt) < 3:
        notes.append("too few replicates for per-region tests (need ≥ 1 residual degree of freedom)")
    comp, amp = amplitude_matrices(a)
    offset, keep = relative_offset(amp["log2_signal"], ref, alt)
    amp["log2_signal_relative"] = amp["log2_signal"] - offset
    if keep.all():
        notes.append("relative offsets: half-sample mode over all regions (iterative refinement not applied)")
    else:
        notes.append(f"relative offsets estimated on {int(keep.sum())}/{len(keep)} presumed-unchanged regions")
    tests = {name: _finish(moderated_t(amp[name], ref, alt), comp, fdr)
             for name in ("log2_auc", "log2_signal", "log2_signal_relative", "log2_enrichment")}
    # Safety net: site-specific changes that trend with signal strength suggest a non-linear
    # normalisation / efficiency difference rather than biology.
    rel = tests["log2_signal_relative"]
    rho = float(stats.spearmanr(rel["estimate"], amp["log2_signal"].mean(axis=0))[0])
    trend = rho
    if abs(rho) > 0.2 and np.isin(rel["call"], ("Higher", "Lower")).mean() > 0.05:
        notes.append(f"relative changes correlate with signal strength (Spearman {rho:+.2f}): part of the site-specific "
                     "calls may reflect non-linear normalisation / IP-efficiency differences")

    # Shape: regions with positive pooled background-corrected mass (condition-blind), optionally
    # restricted to the strongest ones. Statistics use the condition-blind support of each region.
    mass = a.core_mass[:, comp]
    pooled = mass.mean(axis=0)
    keep = (pooled > 0) & (pooled >= np.quantile(pooled, shape_quantile))
    sidx = comp[keep]
    y = a.core_signal[:, sidx]
    e = a.core_edges
    mask = (support_mask(y, KERNEL_SUPPORT_FRACTION) if KERNEL_SUPPORT_FRACTION is not None
            else np.ones(y.shape[1:], bool))
    n_assign = _n_assignments(len(ref), len(alt))
    # Optional condition-blind shape classes on the pooled profiles (descriptive; not used by any test).
    classes = None
    if getattr(c, "shape_classes", False) and len(sidx) >= 60:
        from .shapeclass import shape_classes
        yp = np.clip(np.nan_to_num(y).mean(0), 0, None)
        P = yp / np.maximum(yp.sum(-1, keepdims=True), 1e-12)
        Q = max(1, int(round(c.align_max_shift / float(np.median(np.diff(e))))))
        classes = shape_classes(P, Q, seed=c.seed)
        notes.append(f"shape classes (condition-blind k-mean alignment, descriptive): k = {classes.k}, "
                     f"simplified silhouette {classes.silhouette:.2f}")
    # Processing batches from the manifest: permute labels within batches.
    blocks = None
    batches = [str(a.profiles.samples[i].get("batch", "") or "") for i in np.r_[a.ref_idx, a.alt_idx]]
    if getattr(c, "batch_aware", True) and all(batches) and len(set(batches)) > 1:
        b = np.array(batches)
        n_block = len(label_assignments(len(ref), len(ref) + len(alt), b)[0])
        if n_block < 2:
            blocks = b  # only the observed assignment exists: the tests return NaN (not_tested)
            notes.append("condition and batch coincide (no label assignment within batches exists): per-region shape "
                         "tests are not performed; condition differences cannot be separated from batch")
        elif getattr(c, "batch_permutation", "free") == "within":
            blocks = b
            n_assign = n_block
            notes.append(f"labels permuted within batches ({len(set(batches))} batches; {n_block} assignments)")
        else:
            notes.append(f"batches recorded ({len(set(batches))}); design not confounded; labels permuted freely")
    # With batches that coincide with the condition the tests return NaN: regions are reported as not_tested.
    testable = len(sidx) >= 20 and (n_assign >= 2 or blocks is not None)
    if testable:
        p_kernel = kernel_shape_test(y, mask, ref, alt, KERNEL_SIGMA_BP, KERNEL_S0_QUANTILE, M.centers(e), blocks=blocks)
        tests["shape_kernel"] = _finish(dict(estimate=np.full(len(sidx), np.nan), p=p_kernel, d0=np.nan),
                                        sidx, fdr, omnibus=True)
        notes.append(f"shape_kernel: null pooled over {len(sidx)} regions × {n_assign - 1} non-observed label assignments"
                     + (" (few assignments: conservative)" if n_assign < 10 else ""))
    core = M.core_mask(a.profiles.edges, c.core_flank)
    cov = np.nan_to_num(a.values[:, sidx][..., core])
    # Template fit: effect size of `shift`. Its width scale is not tested: a moderated t on it was
    # miscalibrated under translations (null FWER 0.15–0.50 in the round-2 validation).
    fit = template_fit(y, cov, e, smooth=SHIFT_FIT_SMOOTH, refine=SHIFT_FIT_REFINE)
    if a.profiles.oriented and testable:
        # Calls from the projection score test (pooled permutation null); effect size from the template fit.
        score = shift_score_test(y, ref, alt, M.centers(e), smooth=TEMPLATE_SMOOTH, blocks=blocks)
        est = fit["shift"][alt].mean(0) - fit["shift"][ref].mean(0)
        res = _finish(dict(estimate=score["sign"], p=score["p"], d0=np.nan), sidx, fdr)
        res["estimate"] = est
        res["direction"] = score["sign"]
        tests["shift"] = res
    class_rows = []
    if classes is not None:
        # Summaries of the per-region calls within each class. A class-wide average of all regions'
        # shifts is not reported: most regions do not move, so it is tiny yet "significant" (seen on
        # real yeast MNase data); the direction balance and size among called regions are informative.
        for j in range(classes.k):
            m = classes.labels == j
            row = dict(shape_class=j + 1, regions=int(m.sum()))
            for name in ("shape_kernel", "shift"):
                if name in tests:
                    row[f"{name}_called"] = float(np.isin(tests[name]["call"][m], ("Higher", "Lower", "Changed")).mean())
            if "shift" in tests:
                calls, est = tests["shift"]["call"][m], np.asarray(tests["shift"]["estimate"])[m]
                hi, lo = calls == "Higher", calls == "Lower"
                row.update(shift_higher=int(hi.sum()), shift_lower=int(lo.sum()),
                           shift_higher_share=float(hi.sum() / max(hi.sum() + lo.sum(), 1)),
                           median_abs_shift_called_bp=float(np.median(np.abs(est[hi | lo]))) if (hi | lo).any() else float("nan"))
            class_rows.append(row)
    # Full-window moments on regions with positive mass in every sample (moments are ratios). The
    # support-restricted versions were better on parametric shapes but miscalibrated on real noise
    # (benchmarks/profileseg_region_tuning.py), so they are not used for testing.
    strong = np.all(a.core_mass[:, sidx] > 0, axis=0)
    p_full = M.normalise(a.core_signal[:, sidx[strong]] * M.widths(e))
    shape = {"mad_width": M.mean_abs_deviation(p_full, e)}
    if a.profiles.oriented:
        shape["centroid"] = M.centroid(p_full, e)
    for name, matrix in shape.items():
        tests[name] = _finish(moderated_t(matrix, ref, alt), sidx[strong], fdr)

    # Robust global shape: per replicate, the *median* region-level value; immune to re-weighting.
    # sd_width here (not MAD) so it can be compared with the aggregate sd_width; the median is robust.
    robust = []
    per_region = {"sd_width": M.spread_sd(p_full, e)}
    if a.profiles.oriented:
        per_region["centroid"] = shape["centroid"]
    for name, matrix in per_region.items():
        summary = np.median(matrix, axis=1)
        diff, lo, hi, pw, df = _welch(summary[ref], summary[alt], c.confidence)
        robust.append(dict(metric=f"median_region_{name}", estimate=diff, welch_low=lo, welch_high=hi,
                           welch_p=pw, welch_df=df, n_regions=int(len(sidx))))
    depth = {}
    frags = getattr(a.profiles, "fragments", None)
    if frags is not None and len(sidx):
        use = np.r_[a.ref_idx, a.alt_idx]
        f = np.asarray(frags)[use][:, sidx]
        exact = a.profiles.provenance.get("source") == "bam"
        med = np.nanmedian(f, axis=1)
        low = float(np.mean(np.nanmin(f, axis=0) < LOW_DEPTH_FRAGMENTS))
        names_d = [a.profiles.samples[i]["sample_id"] for i in use]
        depth = dict(source="BAM (exact)" if exact else "bigWig (estimated from the count quantum)",
                     median_fragments=dict(zip(names_d, med.round(1).tolist(), strict=True)),
                     low_depth_share=low, threshold=LOW_DEPTH_FRAGMENTS)
        if np.nanmedian(med) < LOW_DEPTH_FRAGMENTS:
            notes.append(f"low depth: median {np.nanmedian(med):.0f} fragments per region and sample "
                         f"({depth['source']}); per-region shape tests have little power below ~{LOW_DEPTH_FRAGMENTS} — "
                         "a non-significant region is not evidence of no change")
        elif low > 0.25:
            notes.append(f"{low:.0%} of shape-tested regions have fewer than {LOW_DEPTH_FRAGMENTS} fragments in some sample "
                         f"({depth['source']}); power is low there")
    influence = {}
    if len(sidx) >= 20:
        calls = {m: tests[m]["call"] for m in ("shape_kernel", "shift") if m in tests}
        names = [a.profiles.samples[i]["sample_id"] for i in np.r_[a.ref_idx, a.alt_idx]]
        influence = replicate_influence(y, cov, e, ref, alt, calls, fdr, a.profiles.oriented, names)
        notes.extend(influence["flags"])
    return RegionalResult(regions=comp, tests=tests, robust_global=robust, shape_regions=int(len(sidx)), notes=notes,
                          strength_trend=trend, influence=influence, classes=class_rows,
                          class_templates=None if classes is None else classes.templates,
                          class_labels=None if classes is None else classes.labels, depth=depth)


def _n_assignments(n1: int, n2: int) -> int:
    """Distinct replicate-label assignments (an equal split and its mirror count once)."""
    from math import comb
    total = comb(n1 + n2, n1)
    return total // 2 if n1 == n2 else total


def _finish(res: dict, idx: np.ndarray, fdr: float, omnibus: bool = False) -> dict:
    q = bh(res["p"])
    changed = "Changed" if omnibus else None
    call = np.where(~np.isfinite(q), "not_tested",
                    np.where(q < fdr, changed if omnibus else np.where(res["estimate"] > 0, "Higher", "Lower"),
                             "Not_significant"))
    return dict(res, fdr=q, call=call, region_index=idx)


def summary_rows(r: RegionalResult) -> list[dict]:
    rows = []
    for name, t in r.tests.items():
        calls = t["call"]
        rows.append(dict(metric=name, description=REGION_METRICS[name], tested=int(np.sum(calls != "not_tested")),
                         higher=int(np.sum(calls == "Higher")), lower=int(np.sum(calls == "Lower")),
                         changed=int(np.sum(calls == "Changed")),
                         prior_df=float(t["d0"]), median_estimate=float(np.nanmedian(t["estimate"]))))
    return rows
