"""ProfileSeg v1 analysis: QC, amplitude, background-corrected shape, replicate-level effects.

Inference target is *global*: does a metric, summarised over the fixed region set, differ between
conditions? The unit is the biological replicate.

Shape is measured on each replicate's **aggregate** background-subtracted profile. Per-region
truncation of negative noise is avoided because it re-creates the signal-to-noise confound (a
truncated noise floor looks relatively smaller next to a taller peak, so the peak seems
narrower). Before aggregation each region is divided by its *pooled* mass (mean over all
samples), a scale factor that does not depend on condition. Every region therefore contributes
comparably, and within-region amplitude ratios are preserved.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from . import metrics as M
from .data import Config, Profiles

AMPLITUDE_METRICS = ("log2_auc", "log2_centre_auc", "log2_enrichment")
SHAPE_METRICS = ("centroid", "sd_width", "iqr_width", "asymmetry")
NAIVE_METRICS = ("naive_sd_width", "naive_centroid")
SIGNED = {"centroid", "asymmetry", "naive_centroid"}
UNITS = dict(
    log2_auc="log2 FC",
    log2_centre_auc="log2 FC",
    log2_enrichment="log2 FC",
    centroid="bp",
    sd_width="bp",
    iqr_width="bp",
    asymmetry="Δ fraction",
    group_w1="bp",
    naive_sd_width="bp",
    naive_centroid="bp",
)


@dataclass
class Analysis:
    profiles: Profiles
    config: Config
    ref_idx: np.ndarray
    alt_idx: np.ndarray
    core_edges: np.ndarray
    amplitude: dict
    background: np.ndarray
    included: np.ndarray
    exclusion: dict
    p: np.ndarray  # samples × core bins: aggregate background-corrected distributions
    shape: dict  # metric → per-sample values
    naive: dict
    aggregate_negative: np.ndarray
    effects: list = field(default_factory=list)
    region_effects: dict = field(default_factory=dict)
    group: dict = field(default_factory=dict)
    curves: dict = field(default_factory=dict)
    replicate_summaries: dict = field(default_factory=dict)
    diagnostics: dict = field(default_factory=dict)
    flags: list = field(default_factory=list)
    # Retained for per-region (v2) tests: samples × regions × bins, all regions.
    values: np.ndarray | None = None
    core_signal: np.ndarray | None = None  # background-corrected, core window, signed
    complete: np.ndarray | None = None
    core_mass: np.ndarray | None = None


def _shape_metrics(p: np.ndarray, edges: np.ndarray) -> dict:
    return dict(
        centroid=M.centroid(p, edges),
        sd_width=M.spread_sd(p, edges),
        iqr_width=M.iqr_width(p, edges),
        asymmetry=M.asymmetry(p, edges),
    )


def analyse(profiles: Profiles, config: Config, *, bootstrap: int | None = None) -> Analysis:
    cond = profiles.conditions
    ref_idx = np.flatnonzero(cond == config.reference)
    alt_idx = np.flatnonzero(cond == config.alternative)
    if not len(ref_idx) or not len(alt_idx):
        raise ValueError(
            f"need samples for both '{config.reference}' and '{config.alternative}'; found {sorted(set(cond))}"
        )
    flags = []
    others = sorted(set(cond) - {config.reference, config.alternative})
    if others:
        flags.append(f"conditions ignored in this contrast: {others}")
    if len(ref_idx) < 2 or len(alt_idx) < 2:
        flags.append("fewer than 2 replicates in a condition: results are descriptive only")
    use = np.concatenate([ref_idx, alt_idx])
    values = profiles.values[use]
    ref, alt = np.arange(len(ref_idx)), np.arange(len(ref_idx), len(use))
    edges = profiles.edges
    w = M.widths(edges)

    # ---- amplitude on every region without missing bins (fixed, condition-independent set)
    complete = np.all(np.isfinite(values), axis=(0, 2))
    amp = M.amplitude_table(values, edges, config.center_half, config.background_width)

    # ---- background correction on the core window (signed; no per-region truncation)
    bg = M.background(values, edges, config.background_width, config.background_estimator)
    core = M.core_mask(edges, config.core_flank)
    core_edges = np.concatenate([edges[:-1][core], [edges[1:][core][-1]]])
    wc = w[core]
    y = (values - bg[..., None])[..., core]
    mass = np.nansum(y * wc, axis=-1)  # signed background-corrected mass per sample/region

    # ---- condition-independent inclusion: pooled mass only
    reasons = dict(total=int(len(complete)), missing_bins=int((~complete).sum()))
    pooled = mass.mean(axis=0)
    ok = complete & (pooled > 0)
    reasons["nonpositive_pooled_mass"] = int((complete & (pooled <= 0)).sum())
    threshold = (
        float(np.quantile(pooled[ok], config.min_mass_quantile))
        if ok.any() and config.min_mass_quantile > 0
        else 0.0
    )
    low = ok & (pooled < threshold)
    reasons["low_pooled_mass"] = int(low.sum())
    ok &= ~low
    reasons.update(threshold=threshold, included=int(ok.sum()))
    # Informational only (not used for selection): per-condition share of regions whose own mass ≤ 0.
    for name, idx in ((config.reference, ref), (config.alternative, alt)):
        reasons[f"regions_with_nonpositive_mass_in_{name}"] = int(
            np.any(mass[idx][:, ok] <= 0, axis=0).sum()
        )
    if ok.sum() < 3:
        raise ValueError(
            f"only {int(ok.sum())} regions pass inclusion; check background/core windows"
        )

    # Condition-independent per-region scaling. 'shrunk' (default) divides by pooled mass plus its
    # median: weights stay bounded when pooled mass is near zero (low SNR), which otherwise lets a
    # few noisy regions dominate the aggregate. 'equal' = pooled mass only; 'signal' = no scaling.
    if config.region_weighting == "equal":
        scale = pooled[ok]
    elif config.region_weighting == "signal":
        scale = np.full(int(ok.sum()), float(np.median(pooled[ok])))
    elif config.region_weighting == "shrunk_sample":
        # Per-sample bounded scaling: no labels used (permutation stays exact); uniform amplitude
        # changes cancel; region-specific gains re-weight regions much less than pooled scaling.
        m = np.clip(mass[:, ok], 0, None)
        scale = m + np.median(m, axis=1, keepdims=True)
    else:
        scale = pooled[ok] + float(np.median(pooled[ok]))
    scale = np.broadcast_to(scale, mass[:, ok].shape) if np.ndim(scale) == 1 else scale
    Y = y[:, ok] / scale[:, :, None]
    agg = Y.mean(axis=1) * wc  # samples × bins (signed mass)
    negative = np.clip(-agg, 0, None).sum(axis=1) / np.maximum(np.abs(agg).sum(axis=1), M.EPS)
    if config.negative_policy == "exclude" and np.any(negative > config.max_truncated_fraction):
        raise ValueError(
            "aggregate background-corrected profile is substantially negative; "
            "check the background window"
        )
    if np.any(negative > 0.05):
        flags.append(
            f"aggregate profiles carry up to {negative.max():.1%} negative mass after background "
            "correction; shape moments use the signed profile, JSD uses the clipped one"
        )
    # Signed normalisation: zero-mean noise in low-signal tails averages out instead of being
    # clipped into a positive floor (which would again make taller peaks look narrower).
    p = M.normalise(agg)
    bad = ~np.all(np.isfinite(p), axis=-1)
    if bad.any():
        # e.g. MNase arrays, where the flanks carry as much signal as the core: fall back to the
        # clipped profile for those samples and say that aggregate shape metrics are unreliable.
        clipped = M.normalise(np.clip(agg[bad], 0, None))
        if np.any(~np.isfinite(clipped)):
            raise ValueError("an aggregate background-corrected profile has no positive mass; check the background window")
        p = p.copy()
        p[bad] = clipped
        flags.append(f"{int(bad.sum())} sample(s) have non-positive aggregate background-corrected mass; their aggregate "
                     "profile was clipped at zero. Aggregate shape metrics are unreliable here (background window "
                     "carries as much signal as the core); rely on the per-region tests")
    shape = _shape_metrics(p, core_edges)

    # Naive comparator: no background correction, whole window, same pooled scaling.
    raw = np.where(np.isfinite(values[:, ok]), values[:, ok], 0.0)
    raw_pooled = (raw * w).sum(axis=-1).mean(axis=0)
    naive_agg = (raw / np.maximum(raw_pooled, M.EPS)[None, :, None]).mean(axis=1) * w
    pn = M.normalise(naive_agg)
    naive = dict(naive_sd_width=M.spread_sd(pn, edges), naive_centroid=M.centroid(pn, edges))

    a = Analysis(
        profiles=profiles,
        config=config,
        ref_idx=use[ref],
        alt_idx=use[alt],
        core_edges=core_edges,
        amplitude=amp,
        background=bg,
        included=ok,
        exclusion=reasons,
        p=p,
        shape=shape,
        naive=naive,
        aggregate_negative=negative,
        flags=flags,
    )
    _effects(
        a,
        values,
        complete,
        ref,
        alt,
        Y,
        wc,
        raw / np.maximum(raw_pooled, M.EPS)[None, :, None],
        w,
        bootstrap or config.bootstrap,
    )
    diagnostics(a, values, complete, Y, wc, ref, alt)
    a.values, a.core_signal, a.complete, a.core_mass = values, y, complete, mass
    return a


# ---------------------------------------------------------------- inference helpers


def _welch(x, y, confidence):
    nx, ny = len(x), len(y)
    diff = float(np.mean(y) - np.mean(x))
    if nx < 2 or ny < 2:
        return diff, np.nan, np.nan, np.nan, np.nan
    vx, vy = np.var(x, ddof=1) / nx, np.var(y, ddof=1) / ny
    se = np.sqrt(vx + vy)
    if se <= 0:
        return diff, diff, diff, (0.0 if diff else 1.0), np.inf
    df = (vx + vy) ** 2 / (vx**2 / (nx - 1) + vy**2 / (ny - 1))
    t = stats.t.ppf(0.5 + confidence / 2, df)
    return diff, diff - t * se, diff + t * se, float(2 * stats.t.sf(abs(diff) / se, df)), float(df)


def permutation_assignments(
    n_ref: int, n_alt: int, pairs: list | None = None, limit: int = 20000, rng=None
):
    """Distinct replicate-label assignments (within-pair swaps when paired); row 0 is observed."""
    n = n_ref + n_alt
    observed = np.array([False] * n_ref + [True] * n_alt)
    if pairs:
        rows = []
        for flips in itertools.product([False, True], repeat=len(pairs)):
            lab = observed.copy()
            for flip, (i, j) in zip(flips, pairs, strict=True):
                if flip:
                    lab[i], lab[j] = lab[j], lab[i]
            rows.append(lab)
        return np.array(rows), False
    combos = list(itertools.islice(itertools.combinations(range(n), n_alt), limit + 1))
    sampled = len(combos) > limit
    if sampled:
        rng = rng or np.random.default_rng(0)
        combos = [tuple(range(n_ref, n))] + [
            tuple(sorted(rng.choice(n, n_alt, replace=False))) for _ in range(limit - 1)
        ]
    rows = []
    for combo in combos:
        lab = np.zeros(n, dtype=bool)
        lab[list(combo)] = True
        rows.append(lab)
    rows.sort(key=lambda r: not np.array_equal(r, observed))
    return np.array(rows), sampled


def _perm_p(stat_values: np.ndarray):
    """Two-sided p (row 0 observed) and the smallest attainable p for this design."""
    absolute = np.abs(stat_values)
    tol = 1e-9 * max(1.0, float(np.max(absolute)))
    p_obs = float(np.mean(absolute >= absolute[0] - tol))
    floor = float(min(np.mean(absolute >= v - tol) for v in absolute))
    return p_obs, floor


def exact_permutation(summary: np.ndarray, labels: np.ndarray):
    return _perm_p(np.array([summary[lab].mean() - summary[~lab].mean() for lab in labels]))


def _pairs(a: Analysis, ref, alt):
    samples = a.profiles.samples
    ids_ref = [samples[i].get("pairing_id", "") for i in a.ref_idx]
    ids_alt = [samples[i].get("pairing_id", "") for i in a.alt_idx]
    if not all(ids_ref) or sorted(ids_ref) != sorted(ids_alt) or len(set(ids_ref)) != len(ids_ref):
        return None
    return [(int(ref[ids_ref.index(pid)]), int(alt[ids_alt.index(pid)])) for pid in ids_ref]


def _call(diff, lo, hi, n_ref, n_alt):
    if n_ref < 2 or n_alt < 2:
        return "exploratory (fewer than 2 replicates per condition)"
    if lo > 0 or hi < 0:
        return "higher" if diff > 0 else "lower"
    return "no clear change"


def _effects(a: Analysis, values, complete, ref, alt, Y, wc, raw_scaled, w, n_boot):
    c = a.config
    rng = np.random.default_rng(c.seed)
    oriented = a.profiles.oriented
    labels, sampled = permutation_assignments(len(ref), len(alt), _pairs(a, ref, alt), rng=rng)
    if sampled:
        a.flags.append("permutation assignments were sampled (design too large to enumerate)")
    q = (1 - c.confidence) / 2
    n_inc, n_amp = Y.shape[1], int(complete.sum())
    boot_regions = (
        rng.multinomial(n_inc, np.full(n_inc, 1 / n_inc), size=n_boot).astype(float) / n_inc
    )
    boot_amp = rng.multinomial(n_amp, np.full(n_amp, 1 / n_amp), size=n_boot).astype(float) / n_amp
    boot_ref = rng.integers(0, len(ref), size=(n_boot, len(ref)))
    boot_alt = rng.integers(0, len(alt), size=(n_boot, len(alt)))

    def resample_diff(per_sample_boot):  # (samples, boot) → boot differences alt − ref
        r = np.take_along_axis(per_sample_boot[ref].T, boot_ref, axis=1).mean(axis=1)
        al = np.take_along_axis(per_sample_boot[alt].T, boot_alt, axis=1).mean(axis=1)
        return al - r

    # Bootstrap aggregates: resample regions (shared across samples, keeps pairing by region).
    agg_boot = np.einsum("br,srk->sbk", boot_regions, Y) * wc  # samples × boot × bins
    p_boot = M.normalise(agg_boot)
    shape_boot = _shape_metrics(p_boot, a.core_edges)
    naive_boot_agg = np.einsum("br,srk->sbk", boot_regions, raw_scaled) * w
    pn_boot = M.normalise(naive_boot_agg)
    naive_boot = dict(
        naive_sd_width=M.spread_sd(pn_boot, a.profiles.edges),
        naive_centroid=M.centroid(pn_boot, a.profiles.edges),
    )

    summaries, rows = {}, []
    for metric in (*AMPLITUDE_METRICS, *SHAPE_METRICS, *NAIVE_METRICS):
        if metric in SIGNED and not oriented:
            continue
        if metric == "log2_enrichment":
            # centre signal relative to the local background level (both averaged over regions):
            # separates a change of the peak from a global rescaling of signal and background.
            num = a.amplitude["centre_auc"][:, complete]
            den = a.background[:, complete] * (2 * c.center_half)
            summary = np.log2(np.maximum(num.mean(axis=1), M.EPS) / np.maximum(den.mean(axis=1), M.EPS))
            boot = (np.log2(np.maximum(boot_amp @ num.T, M.EPS)) - np.log2(np.maximum(boot_amp @ den.T, M.EPS))).T
        elif metric in AMPLITUDE_METRICS:
            key = "auc" if metric == "log2_auc" else "centre_auc"
            matrix = a.amplitude[key][:, complete]
            summary = np.log2(np.maximum(matrix.mean(axis=1), M.EPS))
            boot = np.log2(np.maximum(boot_amp @ matrix.T, M.EPS)).T
        elif metric in SHAPE_METRICS:
            summary, boot = a.shape[metric], shape_boot[metric]
        else:
            summary, boot = a.naive[metric], naive_boot[metric]
        summaries[metric] = summary
        diff, lo, hi, p_welch, df = _welch(summary[ref], summary[alt], c.confidence)
        bd = resample_diff(boot)
        p_perm, floor = exact_permutation(summary, labels)
        family = (
            "amplitude"
            if metric in AMPLITUDE_METRICS
            else "naive"
            if metric in NAIVE_METRICS
            else "shape"
        )
        rows.append(
            dict(
                metric=metric,
                family=family,
                unit=UNITS[metric],
                estimate=float(diff),
                welch_low=float(lo),
                welch_high=float(hi),
                welch_p=float(p_welch),
                welch_df=float(df),
                boot_low=float(np.quantile(bd, q)),
                boot_high=float(np.quantile(bd, 1 - q)),
                perm_p=p_perm,
                perm_floor=floor,
                n_assignments=len(labels),
                n_ref=len(ref),
                n_alt=len(alt),
                call=_call(diff, lo, hi, len(ref), len(alt)),
            )
        )
    a.replicate_summaries = summaries

    # ---- group-mean distributions and W1 (a distance ≥ 0: permutation, not a CI around 0)
    e = a.core_edges
    p_ref, p_alt = M.normalise(a.p[ref].mean(axis=0)), M.normalise(a.p[alt].mean(axis=0))
    a.group = dict(
        w1=float(M.wasserstein1(p_ref, p_alt, e)),
        w2=float(M.wasserstein2(p_ref, p_alt, e)),
        jsd_bits=float(
            M.jensen_shannon(
                M.normalise(np.clip(p_ref, 0, None)), M.normalise(np.clip(p_alt, 0, None))
            )
        ),
        centroid_shift=float(M.centroid(p_alt, e) - M.centroid(p_ref, e))
        if oriented
        else float("nan"),
        sd_shift=float(M.spread_sd(p_alt, e) - M.spread_sd(p_ref, e)),
    )
    br = M.normalise(
        np.take_along_axis(p_boot[ref].transpose(1, 0, 2), boot_ref[..., None], axis=1).mean(axis=1)
    )
    ba = M.normalise(
        np.take_along_axis(p_boot[alt].transpose(1, 0, 2), boot_alt[..., None], axis=1).mean(axis=1)
    )
    w1_boot = M.wasserstein1(br, ba, e)
    w1_perm = _perm_p(
        np.array(
            [
                M.wasserstein1(
                    M.normalise(a.p[~lab].mean(axis=0)), M.normalise(a.p[lab].mean(axis=0)), e
                )
                for lab in labels
            ]
        )
    )
    rows.append(
        dict(
            metric="group_w1",
            family="shape",
            unit="bp",
            estimate=a.group["w1"],
            welch_low=np.nan,
            welch_high=np.nan,
            welch_p=np.nan,
            welch_df=np.nan,
            boot_low=float(np.quantile(w1_boot, q)),
            boot_high=float(np.quantile(w1_boot, 1 - q)),
            perm_p=w1_perm[0],
            perm_floor=w1_perm[1],
            n_assignments=len(labels),
            n_ref=len(ref),
            n_alt=len(alt),
            call="distance (≥ 0); compare with permutation",
        )
    )
    _assign_roles(a, rows)
    a.effects = rows
    a.curves = dict(
        edges=e,
        p_ref=p_ref,
        p_alt=p_alt,
        diff=p_alt - p_ref,
        diff_low=np.quantile(ba - br, q, axis=0),
        diff_high=np.quantile(ba - br, 1 - q, axis=0),
        raw_edges=a.profiles.edges,
        raw_samples=np.nanmean(values[:, complete], axis=1),
        corrected_samples=a.p,
    )

    # ---- region-level effects (descriptive heterogeneity; per-region profiles are noisy)
    yr = Y[ref].mean(axis=0) * wc  # signed; regions with non-positive mass give NaN
    ya = Y[alt].mean(axis=0) * wc
    pr, pa = M.normalise(yr), M.normalise(ya)
    auc = a.amplitude["auc"][:, a.included]
    pc = 1e-3 * max(float(np.median(auc)), M.EPS)
    nan = np.full(n_inc, np.nan)
    a.region_effects = dict(
        w1=M.wasserstein1(pr, pa, e),
        centroid_shift=(M.centroid(pa, e) - M.centroid(pr, e)) if oriented else nan,
        sd_shift=M.spread_sd(pa, e) - M.spread_sd(pr, e),
        iqr_shift=M.iqr_width(pa, e) - M.iqr_width(pr, e),
        asymmetry_shift=(M.asymmetry(pa, e) - M.asymmetry(pr, e)) if oriented else nan,
        log2fc_auc=np.log2(auc[alt].mean(axis=0) + pc) - np.log2(auc[ref].mean(axis=0) + pc),
    )


# ---------------------------------------------------------------- v1.1: roles, multiplicity, diagnostics


def holm(pvalues: list[float]) -> list[float]:
    """Holm step-down adjusted p-values (NaN entries are left as NaN)."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    idx = np.flatnonzero(np.isfinite(p))
    if not idx.size:
        return out.tolist()
    order = idx[np.argsort(p[idx])]
    m = len(order)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[i]))
        out[i] = running
    return out.tolist()


def primary_metrics(config: Config, oriented: bool) -> tuple[str, str]:
    shape = config.primary_shape
    if shape == "auto":
        shape = "centroid" if oriented else "iqr_width"
    if shape in SIGNED and not oriented:
        raise ValueError(f"primary_shape '{shape}' needs oriented anchors")
    return config.primary_amplitude, shape


def _assign_roles(a: Analysis, rows: list[dict]) -> None:
    """Primary metrics share α by Bonferroni (α/2 each, so the report's family-wise error stays at α);
    secondary metrics are Holm-adjusted together."""
    amp, shape = primary_metrics(a.config, a.profiles.oriented)
    alpha = 1 - a.config.confidence
    n_primary = sum(1 for r in rows if r["metric"] in (amp, shape))
    secondary = [r for r in rows if r["metric"] not in (amp, shape, "group_w1") and r["family"] != "naive"]
    for r, ph in zip(secondary, holm([r["welch_p"] for r in secondary]), strict=True):
        r["p_holm"] = ph
    for r in rows:
        r.setdefault("p_holm", np.nan)
        r.setdefault("p_primary", np.nan)
        if r["metric"] in (amp, shape):
            r["role"] = "primary"
            r["p_primary"] = min(1.0, r["welch_p"] * max(n_primary, 1)) if np.isfinite(r["welch_p"]) else np.nan
            if r["call"] in ("higher", "lower") and not (r["p_primary"] < alpha):
                r["call"] = "no clear change (Bonferroni across primary metrics)"
        elif r["family"] == "naive":
            r["role"] = "comparator"
            r["call"] = "comparator only (no background correction)"
        elif r["metric"] == "group_w1":
            r["role"] = "secondary"
        else:
            r["role"] = "secondary"
            if r["call"] in ("higher", "lower") and not (r["p_holm"] < alpha):
                r["call"] = "no clear change after Holm"
        small = min(r["n_ref"], r["n_alt"]) <= 3
        r["boot_note"] = "descriptive (≤3 replicates per condition)" if small else ""


def diagnostics(a: Analysis, values: np.ndarray, complete: np.ndarray, Y: np.ndarray, wc: np.ndarray,
                ref: np.ndarray, alt: np.ndarray) -> None:
    """Replicate concordance, background-level difference, heterogeneity of region-level shifts."""
    from scipy.stats import rankdata

    d: dict = {}
    names = [a.profiles.samples[i]["sample_id"] for i in np.concatenate([a.ref_idx, a.alt_idx])]
    cond = np.array(["ref"] * len(ref) + ["alt"] * len(alt))
    # 1) replicate concordance: Spearman on per-region log AUC after removing the region effect
    #    (region strength is shared by all samples and would otherwise dominate the correlation).
    auc = np.log2(np.maximum(a.amplitude["auc"][:, complete], M.EPS))
    if auc.shape[1] >= 3 and len(names) > 2:
        resid = auc - auc.mean(axis=0, keepdims=True)
        # ranks + corrcoef always give a samples × samples matrix (spearmanr collapses to a scalar
        # for degenerate input, e.g. noise-free profiles).
        with np.errstate(invalid="ignore", divide="ignore"):
            rho = np.corrcoef(rankdata(resid, axis=1))
        d["concordance"] = dict(samples=names, rho_residual=rho.tolist())
        discordant = []
        for i in range(len(names)):
            same = [rho[i, j] for j in range(len(names)) if j != i and cond[j] == cond[i]]
            other = [rho[i, j] for j in range(len(names)) if cond[j] != cond[i]]
            if same and other and np.mean(same) < np.mean(other):
                discordant.append(names[i])
        # Informational: residual correlation reflects region-specific variation only (uniform
        # amplitude changes and pure shape changes leave it untouched), so it is reported in QC
        # rather than used to judge effects.
        d["discordant_replicates"] = discordant
    # Direction consistency of each called effect across individual replicates.
    for r in a.effects:
        summary = a.replicate_summaries.get(r["metric"])
        if summary is None or not np.isfinite(r["estimate"]):
            r["replicate_consistency"] = ""
            continue
        mid = (summary[ref].mean() + summary[alt].mean()) / 2
        sign = np.sign(r["estimate"])
        ok_alt = int(np.sum(np.sign(summary[alt] - mid) == sign))
        ok_ref = int(np.sum(np.sign(summary[ref] - mid) == -sign))
        r["replicate_consistency"] = f"{ok_ref}/{len(ref)} ref, {ok_alt}/{len(alt)} alt on expected side"
        if r["role"] == "primary" and r["call"] in ("higher", "lower") and (ok_ref < len(ref) or ok_alt < len(alt)):
            a.flags.append(
                f"{r['metric']} is called {r['call']} but not every replicate lies on the expected side "
                f"({r['replicate_consistency']}); the effect may be driven by individual replicates"
            )
    # 2) background level difference between conditions
    bg = a.background[:, complete]
    ratio = float(np.log2(np.maximum(bg[alt].mean(), M.EPS) / np.maximum(bg[ref].mean(), M.EPS)))
    d["background_log2_ratio"] = ratio
    if abs(ratio) > 0.25:
        a.flags.append(
            f"local background differs between conditions (log2 ratio {ratio:+.2f}): amplitude contrasts "
            "may reflect normalisation or IP efficiency; compare log2_enrichment, or use spike-in/count models"
        )
    # 3) heterogeneity: region-level shifts versus a within-condition replicate-split noise floor
    e = a.core_edges

    def region_shifts(ga, gb):
        pa_ = M.normalise(Y[ga].mean(axis=0) * wc)
        pb_ = M.normalise(Y[gb].mean(axis=0) * wc)
        out = {"sd": M.spread_sd(pb_, e) - M.spread_sd(pa_, e)}
        if a.profiles.oriented:
            out["centroid"] = M.centroid(pb_, e) - M.centroid(pa_, e)
        return out

    between = region_shifts(ref, alt)
    splits = [g for g in (ref, alt) if len(g) >= 2]
    if splits:
        within = [region_shifts(g[: len(g) // 2], g[len(g) // 2:]) for g in splits]
        het = {}
        for key in between:
            noise = float(np.nanmedian(np.concatenate([np.abs(w[key]) for w in within])))
            signal = float(np.nanmedian(np.abs(between[key])))
            het[key] = dict(median_abs_between=signal, median_abs_within_split=noise,
                            ratio=signal / noise if noise > 0 else np.nan)
        d["heterogeneity"] = het
        calls = {r["metric"]: r for r in a.effects}
        width_changed = any(calls.get(m, {}).get("call") in ("higher", "lower") for m in ("sd_width", "iqr_width"))
        centroid_flat = calls.get("centroid", {}).get("call") not in ("higher", "lower")
        if "centroid" in het and width_changed and centroid_flat and het["centroid"]["ratio"] > 1.5:
            a.flags.append(
                f"aggregate width change with large region-level shifts in both directions "
                f"(median |shift| {het['centroid']['median_abs_between']:.0f} bp vs within-condition "
                f"{het['centroid']['median_abs_within_split']:.0f} bp): possibly bidirectional shifts, "
                "not broadening; see per-region tests"
            )
    # 4) re-weighting caveat: region-specific amplitude changes re-weight regions in the aggregate,
    #    which can move aggregate shape when the affected regions have a different shape.
    lfc = np.log2(np.maximum(a.amplitude["auc"][:, complete], M.EPS))
    between_lfc = lfc[alt].mean(axis=0) - lfc[ref].mean(axis=0)
    spread_between = float(np.std(between_lfc - np.median(between_lfc)))
    split_sd = []
    for g in (ref, alt):
        if len(g) >= 2:
            w_lfc = lfc[g[len(g) // 2:]].mean(axis=0) - lfc[g[: len(g) // 2]].mean(axis=0)
            split_sd.append(float(np.std(w_lfc - np.median(w_lfc))))
    if split_sd:
        ratio = spread_between / max(np.mean(split_sd), M.EPS)
        d["amplitude_heterogeneity_ratio"] = ratio
        shape_called = any(r["call"] in ("higher", "lower") for r in a.effects if r["family"] == "shape")
        if ratio > 1.5 and shape_called:
            a.flags.append(
                f"region-specific amplitude changes are {ratio:.1f}× the within-condition noise; the aggregate "
                "shape change may partly reflect re-weighting of regions rather than a change of shape — "
                "check per-region shape tests"
            )
    a.diagnostics = d
