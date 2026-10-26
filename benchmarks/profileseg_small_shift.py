"""Why are 50 bp shifts hard per region, and what helps? Diagnostics and candidate statistics.

All candidates keep replicates as the unit and use the same pooled permutation null as
``kernel_shape_test`` (null pooled over regions and non-observed label assignments, label-invariant s0).

  kernel_s{30,60,120}  current omnibus kernel statistic at three scales
  kernel_multi         max over the three scales of rank-normalised statistics (same permutations)
  proj_shift           1-df projection of each sample's profile on the template derivative
                       (locally most powerful direction for a small shift)
  proj_2df             projections on the derivative and on (x − μ)·derivative (shift + width)
  *_ihw                the same p-values with weighted BH: weights from a label-invariant
                       sharpness × depth score (independent hypothesis weighting, grouped)

Diagnostic: per-sample shift estimates (projection) → replicate SD per region, and the implied
z-score for a 50 / 100 bp shift with n replicates per group.

    python benchmarks/profileseg_small_shift.py --out DIR --templates T.npz --semireal CACHE.npz
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

import profileseg_region_methods as RM
import profileseg_region_tuning as T
from genomecanvas.profileseg.regional import bh, moderated_t
from genomecanvas.profileseg.simulate import SimSettings

SOURCES = ("param", "templates", "semireal")
DESIGNS = (3, 2)
SCENARIOS = ("null", "mix_shift50", "mix_shift100", "mix_broaden")
FDR = 0.05


def _splits(n1, n):
    return [s for s in itertools.combinations(range(n), n1) if n1 != n - n1 or 0 in s]


def pooled_p(F: np.ndarray, obs_i: int) -> np.ndarray:
    obs = F[obs_i]
    null = np.sort(np.delete(F, obs_i, axis=0).ravel())
    return (1 + null.size - np.searchsorted(null, obs, side="left")) / (1 + null.size)


def kernel_F(p, n1, centers, sigma, s0q=0.5):
    n = p.shape[0]
    K = np.exp(-0.5 * ((centers[:, None] - centers[None]) / sigma) ** 2)
    dev = p - p.mean(0)
    s0 = np.quantile(np.einsum("srb,bc,src->r", dev, K, dev) / (n - 1), s0q)
    F = []
    splits = _splits(n1, n)
    for s in splits:
        a = np.array(s); b = np.setdiff1d(np.arange(n), a)
        ma, mb = p[a].mean(0), p[b].mean(0)
        D = np.einsum("rb,bc,rc->r", ma - mb, K, ma - mb)
        res = np.concatenate([p[a] - ma, p[b] - mb])
        W = np.einsum("srb,bc,src->r", res, K, res) / (n - 2)
        F.append(D / (W + s0))
    return np.array(F), splits.index(tuple(range(n1)))


def projections(p, centers, smooth=5):
    """Per-sample coefficients on the template derivative (≈ −shift) and on (x−μ)·derivative."""
    Tm = p.mean(0)
    half = (smooth + 1) // 2
    k = np.r_[np.arange(1, half + 1), np.arange(half - 1, 0, -1)].astype(float); k /= k.sum()
    Tm = np.apply_along_axis(lambda r: np.convolve(r, k, mode="same"), 1, Tm)
    h = float(np.median(np.diff(centers)))
    dT = np.gradient(Tm, h, axis=-1)
    mu = (Tm * centers).sum(1, keepdims=True) / np.maximum(Tm.sum(1, keepdims=True), 1e-12)
    g1, g2 = dT, (centers - mu) * dT
    dev = p - Tm
    z1 = (dev * g1).sum(-1) / np.maximum((g1 * g1).sum(-1), 1e-12)  # shift ≈ −z1
    z2 = (dev * g2).sum(-1) / np.maximum((g2 * g2).sum(-1), 1e-12)
    return -z1, z2


def proj_F(z: np.ndarray, n1: int, s0q=0.5):
    """z: samples × regions × k. F = ||Δ mean||² / (within variance + s0), pooled null."""
    n = z.shape[0]
    dev = z - z.mean(0)
    s0 = np.quantile((dev ** 2).sum((0, 2)) / (n - 1), s0q)
    F, splits = [], _splits(n1, n)
    for s in splits:
        a = np.array(s); b = np.setdiff1d(np.arange(n), a)
        d = z[a].mean(0) - z[b].mean(0)
        W = (((z[a] - z[a].mean(0)) ** 2).sum((0, 2)) + ((z[b] - z[b].mean(0)) ** 2).sum((0, 2))) / (n - 2)
        F.append((d ** 2).sum(-1) / (W + s0))
    return np.array(F), splits.index(tuple(range(n1)))


def weighted_bh(p: np.ndarray, score: np.ndarray, groups: int = 5, fdr: float = FDR) -> np.ndarray:
    """Grouped weighted BH (Genovese et al. 2006 / IHW-style): weights ∝ estimated non-null share
    per score group (Storey π0 at λ = 0.5), normalised to mean 1. Score is label-invariant."""
    q = np.quantile(score, np.linspace(0, 1, groups + 1))
    g = np.clip(np.searchsorted(q, score, side="right") - 1, 0, groups - 1)
    w = np.ones_like(p)
    for j in range(groups):
        m = g == j
        if m.sum() < 10:
            continue
        pi0 = min(1.0, np.mean(p[m] > 0.5) / 0.5)
        w[m] = max(1 - pi0, 0.02)
    w = w / w.mean()
    pw = np.minimum(p / np.maximum(w, 1e-6), 1.0)
    return bh(pw) < fdr


def evaluate(d: RM.Data, n: int) -> dict:
    y = np.clip(d.y, 0, None)
    p = y / np.maximum(y.sum(-1, keepdims=True), 1e-12)
    c = d.centers
    res, Fs = {}, {}
    for sig in (30, 60, 120):
        F, o = kernel_F(p, n, c, sig)
        Fs[sig] = (F, o)
        res[f"kernel_s{sig}"] = pooled_p(F, o)
    # multi-scale: rank-normalise each scale over all (assignment, region) values, take the max
    ranks = []
    for sig, (F, o) in Fs.items():
        flat = F.ravel()
        r = np.empty_like(flat); r[np.argsort(flat)] = np.arange(flat.size)
        ranks.append(r.reshape(F.shape) / flat.size)
    res["kernel_multi"] = pooled_p(np.max(ranks, axis=0), Fs[60][1])
    shift, scale = projections(p, c)
    F, o = proj_F(shift[..., None], n)
    res["proj_shift"] = pooled_p(F, o)
    F, o = proj_F(shift[..., None], n, s0q=0.1)
    res["proj_shift_s0low"] = pooled_p(F, o)
    res["proj_shift_modt"] = moderated_t(shift, np.arange(n), np.arange(n, 2 * n))["p"]
    F, o = proj_F(np.stack([shift / np.std(shift), scale / np.std(scale)], -1), n)
    res["proj_2df"] = pooled_p(F, o)
    # label-invariant priority score: depth × sharpness of the pooled template
    Tm = p.mean(0)
    sharp = np.sqrt((np.gradient(Tm, axis=-1) ** 2).sum(-1))
    depth = d.y.sum(-1).mean(0)
    score = np.log(np.maximum(sharp, 1e-12)) + 0.5 * np.log(np.maximum(depth, 1e-12))
    out = {}
    for k, pv in res.items():
        out[k] = RM.score(pv, d.changed)
        called = weighted_bh(pv, score)
        tr = d.changed; nc = int(called.sum())
        out[k + "_ihw"] = dict(called=float(called.mean()), any=float(called.any()),
                               fdp=float((called & ~tr).sum() / nc) if nc else 0.0,
                               power=float((called & tr).sum() / tr.sum()) if tr.any() else float("nan"))
    # diagnostic: replicate SD of per-sample shift estimates in unchanged regions
    unchanged = ~d.changed
    sd = np.sqrt((np.var(shift[:n, unchanged], axis=0, ddof=1) + np.var(shift[n:, unchanged], axis=0, ddof=1)) / 2)
    ch = d.changed
    est = shift[n:].mean(0) - shift[:n].mean(0)
    sd_all = np.sqrt((np.var(shift[:n], axis=0, ddof=1) + np.var(shift[n:], axis=0, ddof=1)) / 2)
    z_real = np.abs(est) / (sd_all * np.sqrt(2 / n))
    out["_diag"] = dict(median_rep_sd_bp=float(np.median(sd)),
                        est_changed=float(np.median(est[ch])) if ch.any() else float("nan"),
                        z_realised_changed=float(np.median(z_real[ch])) if ch.any() else float("nan"),
                        z_realised_unchanged=float(np.median(z_real[~ch])),
                        z50_median=float(np.median(50 / (sd * np.sqrt(2 / n)))),
                        z100_median=float(np.median(100 / (sd * np.sqrt(2 / n)))),
                        share_z50_gt3=float(np.mean(50 / (sd * np.sqrt(2 / n)) > 3)))
    return out


def run(job):
    source, n, scenario, rep, templates_path, cache = job
    seed = 500_000 + 10_000 * SOURCES.index(source) + 1000 * SCENARIOS.index(scenario) + 100 * n + rep
    shift = 50.0 if scenario == "mix_shift50" else 100.0
    base = scenario.replace("50", "").replace("100", "")
    if source == "semireal":
        if cache not in T._SEMI:
            from semireal import SemiReal
            T._SEMI[cache] = SemiReal(Path(cache))
        from genomecanvas.profileseg import metrics as M
        from genomecanvas.profileseg.analysis import analyse
        from genomecanvas.profileseg.data import Config
        change = {"null": "null", "mix_shift": "shift", "mix_broaden": "broaden"}[base]
        prof, changed = T._SEMI[cache].experiment(change, 300, 0.2, seed, n_per_group=n, shift_bp=shift)
        config = Config(reference="WT", alternative="KO")
        a = analyse(prof, config, bootstrap=0)
        ok = a.complete & (np.nanmean(a.core_mass, axis=0) > 0)
        d = RM.Data(np.nan_to_num(a.core_signal)[:, ok], None, None, M.centers(a.core_edges), changed[ok])
    else:
        templates = np.load(templates_path)["templates"] if source == "templates" else None
        d = RM.prepare(base, SimSettings(n_regions=300, n_ref=n, n_alt=n, shift_bp=shift), templates, seed)
    return source, n, scenario, evaluate(d, n)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--templates", required=True)
    ap.add_argument("--semireal", required=True)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    jobs = [(s, n, sc, r, args.templates, args.semireal) for s in SOURCES for n in DESIGNS for sc in SCENARIOS
            for r in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for s, n, sc, res in pool.map(run, jobs):
            for m, v in res.items():
                agg.setdefault(f"{s}:{n}v{n}", {}).setdefault(sc, {}).setdefault(m, []).append(v)
    table = {cell: {sc: {m: {k: float(np.nanmean([r[k] for r in rows])) for k in rows[0]} for m, rows in ms.items()}
                    for sc, ms in d.items()} for cell, d in agg.items()}
    (out / "small_shift.json").write_text(json.dumps(table, indent=1))
    cells = sorted(table)
    methods = [m for m in next(iter(table.values()))["null"] if not m.startswith("_")]
    lines = ["# Small shifts: candidate statistics", "", "Power at 50 bp (20 % of regions shifted), BH 5 %; "
             "worst null FWER, worst FDP over all mixtures, and power at 100 bp / broaden for reference.", "",
             "| method | " + " | ".join(cells) + " | worst null FWER | worst FDP | mean power 100 bp | mean power broaden |",
             "|---|" + "---|" * (len(cells) + 4)]
    for m in methods:
        p50 = " | ".join(f"{table[c]['mix_shift50'][m]['power']:.2f}" for c in cells)
        fw = max(table[c]["null"][m]["any"] for c in cells)
        fdp = max(table[c][sc][m]["fdp"] for c in cells for sc in ("mix_shift50", "mix_shift100", "mix_broaden"))
        p100 = np.mean([table[c]["mix_shift100"][m]["power"] for c in cells])
        pb = np.mean([table[c]["mix_broaden"][m]["power"] for c in cells])
        lines.append(f"| {m} | {p50} | {fw:.2f} | {fdp:.3f} | {p100:.2f} | {pb:.2f} |")
    lines += ["", "## Information: replicate SD of per-sample shift estimates (unchanged regions)", "",
              "| cell | median replicate SD (bp) | nominal z, 50 bp | estimated shift, 50 bp truth (bp) | realised z, shifted regions | realised z, unshifted |",
              "|---|---|---|---|---|---|"]
    for c in cells:
        dg = table[c]["mix_shift50"]["_diag"]
        lines.append(f"| {c} | {dg['median_rep_sd_bp']:.0f} | {dg['z50_median']:.2f} | {dg['est_changed']:.1f} | "
                     f"{dg['z_realised_changed']:.2f} | {dg['z_realised_unchanged']:.2f} |")
    (out / "SMALL_SHIFT.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
