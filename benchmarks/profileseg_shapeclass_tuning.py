"""Shape classes and anchor alignment (FunChIP ideas): do they help, and do they keep calibration?

Tuning seeds (1 300 000+). Variants, end to end through ``analyse`` + ``per_region_tests``:

  V0  current v2.2: global pooled nulls, no alignment
  V1  pooled nulls stratified by condition-blind shape classes (rejected; reproduced here with the
      low-level ``strata`` arguments)
  V2  anchors re-centred on condition-blind class templates (no strata)
  V3  strata + alignment
  plus `kernel_d0d1` (V0 data): kernel statistic on profiles and on their first derivatives,
  combined by the maximum of rank-normalised statistics (same permutations)

Scenarios: null, other no-shape changes, 20 % mixtures (shift 100 / 50 bp, broaden), every region
shifted by 50 bp (aggregate detection), each "+jitter" variant adding a region-level anchor error
(uniform ±100 bp, identical for all samples of the region).

    python benchmarks/profileseg_shapeclass_tuning.py --out DIR --templates T.npz --semireal CACHE.npz
"""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, Profiles
from genomecanvas.profileseg.regional import (KERNEL_S0_QUANTILE, KERNEL_SIGMA_BP, TEMPLATE_SMOOTH, bh,
                                               kernel_shape_test, per_region_tests, pooled_permutation_p,
                                               shift_score_test)
from genomecanvas.profileseg.shapeclass import align_profiles, shape_classes
from genomecanvas.profileseg.simulate import SimSettings, simulate
from profileseg_region_validation import _SEMI, make
from profileseg_small_shift import kernel_F

SOURCES = ("param", "templates", "semireal")
DESIGNS = (3, 2)
SCENARIOS = ("null", "null+jitter", "other_amplitude", "mix_shift", "mix_shift+jitter", "mix_shift50",
             "mix_broaden", "shift_all50", "shift_all50+jitter")
VARIANTS = {"V0": (False, False), "V1": (True, False), "V2": (False, True), "V3": (True, True)}
SEED_BASE = 1_300_000


def jitter(prof: Profiles, rng, max_bp=100.0) -> Profiles:
    w = float(np.median(np.diff(prof.edges)))
    q = rng.integers(-int(max_bp / w), int(max_bp / w) + 1, len(prof.regions))
    v = prof.values.copy()
    for r, s in enumerate(q):
        if s:
            row = v[:, r]
            out = np.empty_like(row)
            if s > 0:
                out[:, s:], out[:, :s] = row[:, :-s], row[:, :1]
            else:
                out[:, :s], out[:, s:] = row[:, -s:], row[:, -1:]
            v[:, r] = out
    return Profiles(v, prof.edges, prof.samples, prof.regions, prof.oriented)


def build(source, n, scenario, seed, tpl, cache):
    base = scenario.replace("+jitter", "")
    if base == "shift_all50":
        if source == "semireal":
            from semireal import SemiReal
            if cache not in _SEMI:
                _SEMI[cache] = SemiReal(Path(cache))
            prof, changed = _SEMI[cache].experiment("shift", 300, 1.0, seed, n_per_group=n, shift_bp=50.0)
        else:
            templates = np.load(tpl)["templates"] if source == "templates" else None
            prof, _ = simulate("shift", SimSettings(n_regions=300, n_ref=n, n_alt=n, shift_bp=50.0), seed=seed,
                               templates=templates)
            changed = np.ones(300, bool)
    else:
        prof, changed, _ = make(source, n, base, seed, tpl, cache)
    if scenario.endswith("+jitter"):
        prof = jitter(prof, np.random.default_rng(seed + 1))
    return prof, changed


def score(calls, truth):
    nc = int(calls.sum())
    return dict(called=float(calls.mean()), any=float(calls.any()),
                fdp=float((calls & ~truth).sum() / nc) if nc else 0.0,
                power=float((calls & truth).sum() / truth.sum()) if truth.any() else float("nan"))


def run(job):
    source, n, scenario, rep, tpl, cache = job
    seed = SEED_BASE + 10_000 * SOURCES.index(source) + 1000 * SCENARIOS.index(scenario) + 100 * n + rep
    prof, changed = build(source, n, scenario, seed, tpl, cache)
    out = {}
    for name, (strata, align) in VARIANTS.items():
        cfg = Config(reference="WT", alternative="KO", align_anchors=align)
        p = prof
        if align:
            p, _, _ = align_profiles(prof, cfg.background_width, cfg.core_flank, cfg.align_max_shift, seed=cfg.seed)
        a = analyse(p, cfg, bootstrap=0)
        reg = per_region_tests(a)
        res = {}
        if strata:  # stratified pooled nulls (option removed from per_region_tests after this experiment)
            idx = reg.tests["shape_kernel"]["region_index"]
            y = a.core_signal[:, idx]
            c = M.centers(a.core_edges)
            yp = np.clip(np.nan_to_num(y).mean(0), 0, None)
            P = yp / np.maximum(yp.sum(-1, keepdims=True), 1e-12)
            sc = shape_classes(P, int(round(cfg.align_max_shift / float(np.median(np.diff(a.core_edges))))), seed=cfg.seed)
            lab = sc.labels if sc.k > 1 else None
            ref, alt = np.arange(n), np.arange(n, 2 * n)
            pk = kernel_shape_test(y, np.ones(y.shape[1:], bool), ref, alt, KERNEL_SIGMA_BP, KERNEL_S0_QUANTILE, c, lab)
            ps = shift_score_test(y, ref, alt, c, smooth=TEMPLATE_SMOOTH, strata=lab)["p"]
            for m, pv in (("shape_kernel", pk), ("shift", ps)):
                res[m] = score(bh(pv) < 0.05, changed[idx])
            res["k"] = sc.k
        else:
            for m in ("shape_kernel", "shift"):
                t = reg.tests[m]
                calls = np.isin(t["call"], ("Higher", "Lower", "Changed"))
                res[m] = score(calls, changed[t["region_index"]])
        eff = {e["metric"]: e for e in a.effects}
        res["aggregate_centroid"] = float(eff.get("centroid", {}).get("call") in ("higher", "lower"))
        res.setdefault("k", 1)
        out[name] = res
        if name == "V0":  # derivative-augmented kernel on the same data
            y = np.clip(np.nan_to_num(a.core_signal[:, reg.tests["shape_kernel"]["region_index"]]), 0, None)
            pn = y / np.maximum(y.sum(-1, keepdims=True), 1e-12)
            c = M.centers(a.core_edges)
            F0, o = kernel_F(pn, n, c, 60.0)
            F1, _ = kernel_F(np.gradient(pn, axis=-1), n, c, 60.0)
            ranks = []
            for F in (F0, F1):
                flat = F.ravel(); r = np.empty_like(flat); r[np.argsort(flat)] = np.arange(flat.size)
                ranks.append(r.reshape(F.shape) / flat.size)
            pv = pooled_permutation_p(np.max(ranks, axis=0), o)
            out["kernel_d0d1"] = {"shape_kernel": score(bh(pv) < 0.05, changed[reg.tests["shape_kernel"]["region_index"]])}
    return source, n, scenario, out


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
            for v, r in res.items():
                agg.setdefault(v, {}).setdefault(f"{s}:{n}v{n}", {}).setdefault(sc, []).append(r)
    mean = lambda rows, m, k: float(np.nanmean([r[m][k] if isinstance(r[m], dict) else r[m] for r in rows]))
    table = {}
    lines = ["# Shape classes and anchor alignment (tuning seeds)", ""]
    no_shape = ("null", "null+jitter", "other_amplitude")
    mixes = ("mix_shift", "mix_shift+jitter", "mix_shift50", "mix_broaden")
    for v in agg:
        tests = ("shape_kernel",) if v == "kernel_d0d1" else ("shape_kernel", "shift")
        for m in tests:
            cells = sorted(agg[v])
            fw = max(mean(agg[v][c]["null"], m, "any") for c in cells)
            fs = max(mean(agg[v][c][sc], m, "called") for c in cells for sc in no_shape)
            fdp = max(mean(agg[v][c][sc], m, "fdp") for c in cells for sc in mixes)
            pw = {sc: np.mean([mean(agg[v][c][sc], m, "power") for c in cells]) for sc in mixes}
            table[f"{v}:{m}"] = dict(max_null_fwer=fw, max_false_share=fs, max_fdp=fdp, power=pw)
    lines += ["## Per region (worst cell for errors; mean over cells for power)", "",
              "| variant · test | max null FWER | max false share | max FDP | power shift100 | shift100 + jitter | shift50 | broaden |",
              "|---|---|---|---|---|---|---|---|"]
    for k, t in table.items():
        p = t["power"]
        lines.append(f"| {k} | {t['max_null_fwer']:.2f} | {t['max_false_share']:.4f} | {t['max_fdp']:.3f} | "
                     f"{p['mix_shift']:.2f} | {p['mix_shift+jitter']:.2f} | {p['mix_shift50']:.2f} | {p['mix_broaden']:.2f} |")
    lines += ["", "## Aggregate centroid call (share of runs)", "",
              "| variant | " + " | ".join(sorted(agg["V0"])) + " |", "|---|" + "---|" * len(agg["V0"])]
    for v in VARIANTS:
        for sc in ("null", "shift_all50", "shift_all50+jitter"):
            cells = sorted(agg[v])
            lines.append(f"| {v} {sc} | " + " | ".join(f"{mean(agg[v][c][sc], 'aggregate_centroid', 'x') if False else np.mean([r['aggregate_centroid'] for r in agg[v][c][sc]]):.2f}" for c in cells) + " |")
    lines += ["", "## Number of shape classes chosen (V1, mean)", "",
              " | ".join(f"{c}: {np.mean([r['k'] for sc in SCENARIOS for r in agg['V1'][c][sc]]):.1f}" for c in sorted(agg['V1']))]
    (out / "shapeclass_tuning.json").write_text(json.dumps(table, indent=1, default=float))
    (out / "SHAPECLASS_TUNING.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
