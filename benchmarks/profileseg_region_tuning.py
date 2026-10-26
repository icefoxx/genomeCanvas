"""Tune the per-region shape methods on tuning seeds (validation uses different seeds).

Sources: parametric simulation, real-shape templates, and semi-real data (real replicate noise,
`semireal.py`); designs 3 vs 3 and 2 vs 2. Scenarios: null, amplitude-type changes without shape change, and 20 % mixtures of
shift (100 bp) or broadening (×1.5).

Selection rule (fixed before running):
  E (kernel test): among variants whose complete-null FWER is ≤ 0.10, whose false-call share is
     ≤ 0.005 in every no-shape scenario of every source, and whose mean FDP is ≤ 0.10 in every
     mixture cell (this last condition was added after the first tuning run, before any validation
     run, when the first pick had FDP 0.29 in one cell), take the highest mean power over
     {mix_shift, mix_broaden} × {param, templates, semireal} × {3v3, 2v2}.
  A (support moments): same constraint; centroid power on mix_shift + MAD power on mix_broaden.
  C2 (template fit): same constraint on the shift test; smallest mean |bias| of the shift estimate,
     power as tie-break.

    python benchmarks/profileseg_region_tuning.py --out DIR --templates T.npz --semireal CACHE.npz
"""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

import profileseg_region_methods as RM
from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.regional import (kernel_shape_test, moderated_t, support_mask, support_moments,
                                               template_fit)
from genomecanvas.profileseg.simulate import SimSettings
from semireal import SemiReal

SCENARIOS = {"param": ("null", "amplitude", "site_amplitude", "mix_shift", "mix_broaden"),
             "templates": ("null", "amplitude", "site_amplitude", "mix_shift", "mix_broaden"),
             "semireal": ("null", "amplitude", "global_scale", "mix_shift", "mix_broaden")}
E_GRID = [(sig, frac, q) for sig in (30, 60, 120) for frac in (0.1, 0.2, None) for q in (0.1, 0.5, 0.9)]
A_GRID = (0.1, 0.2, 0.3)
C2_GRID = (1, 3, 5)
TRUE_SHIFT = 100.0
DESIGNS = ("3v3", "2v2")
KEYS = [f"{s}:{d}" for s in ("param", "templates", "semireal") for d in DESIGNS]
_SEMI = {}


def semireal_data(cache: str, scenario: str, seed: int, n: int) -> RM.Data:
    if cache not in _SEMI:
        _SEMI[cache] = SemiReal(Path(cache))
    sr = _SEMI[cache]
    change, frac = (scenario[4:], 0.2) if scenario.startswith("mix_") else (scenario, 0.2)
    prof, changed = sr.experiment(change, 300, frac, seed, n_per_group=n)
    config = Config(reference="WT", alternative="KO")
    a = analyse(prof, config, bootstrap=0)
    ok = a.complete & (np.nanmean(a.core_mass, axis=0) > 0)
    core = M.core_mask(prof.edges, config.core_flank)
    cov = np.nan_to_num(prof.values[:, :, core])
    return RM.Data(np.nan_to_num(a.core_signal)[:, ok], cov[:, ok], np.zeros_like(cov[:, ok], dtype=int),
                   M.centers(a.core_edges), changed[ok])


def run(job):
    source, design, scenario, rep, templates_path, cache = job
    n = int(design[0])
    seed = 100_000 + 1000 * list(SCENARIOS[source]).index(scenario) + 100 * n + rep
    if source == "semireal":
        d = semireal_data(cache, scenario, seed, n)
    else:
        templates = np.load(templates_path)["templates"] if source == "templates" else None
        d = RM.prepare(scenario, SimSettings(n_regions=300, n_ref=n, n_alt=n), templates, seed)
    ref, alt = np.arange(n), np.arange(n, 2 * n)
    edges = np.r_[d.centers - 10, d.centers[-1] + 10]
    out = {}
    masks = {frac: support_mask(d.y, frac) for frac in (0.1, 0.2, 0.3)}
    masks[None] = np.ones(d.y.shape[1:], bool)
    for sig, frac, q in E_GRID:
        out[f"E|{sig}|{frac}|{q}"] = (kernel_shape_test(d.y, masks[frac], ref, alt, sig, q, d.centers), None)
    for frac in A_GRID:
        cs, ms = support_moments(d.y, edges, masks[frac])
        out[f"Ac|{frac}"] = (moderated_t(cs, ref, alt)["p"], None)
        out[f"Am|{frac}"] = (moderated_t(ms, ref, alt)["p"], None)
    for sm in C2_GRID:
        tf = template_fit(d.y, d.cov, edges, smooth=sm)
        r = moderated_t(tf["shift"], ref, alt)
        out[f"C2|{sm}"] = (r["p"], r["estimate"])
    res = {}
    for k, (p, est) in out.items():
        s = RM.score(p, d.changed)
        if est is not None and d.changed.any():
            s["bias"] = float(np.median(np.asarray(est)[d.changed]) - TRUE_SHIFT) if "shift" in scenario else np.nan
        res[k] = s
    return f"{source}:{design}", scenario, res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--templates", required=True)
    ap.add_argument("--semireal", required=True)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    jobs = [(src, des, sc, rep, args.templates, args.semireal) for src, scs in SCENARIOS.items() for des in DESIGNS
            for sc in scs for rep in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for src, sc, res in pool.map(run, jobs):
            for k, s in res.items():
                agg.setdefault(k, {}).setdefault(src, {}).setdefault(sc, []).append(s)
    table = {k: {src: {sc: {m: float(np.nanmean([r.get(m, np.nan) for r in rows])) for m in ("called", "any", "power", "fdp", "bias")}
                       for sc, rows in d2.items()} for src, d2 in d1.items()} for k, d1 in agg.items()}
    (out / "tuning.json").write_text(json.dumps(table, indent=1))

    def admissible(k):
        for src, scs in table[k].items():
            for sc, t in scs.items():
                if sc == "null" and t["any"] > 0.10:
                    return False
                if not sc.startswith("mix_") and t["called"] > 0.005:
                    return False
                if sc.startswith("mix_") and t["fdp"] > 0.10:  # amended before validation: matches criterion V3
                    return False
        return True

    def power(k, scen):
        return float(np.mean([table[k][src][sc]["power"] for src in table[k] for sc in scen]))

    lines = ["# Per-region tuning", "", "| variant | admissible | mix_shift power (" + "/".join(KEYS) + ") | mix_broaden power | max null FWER | max false share | mean FDP (mix) | bias (bp) |",
             "|---|---|---|---|---|---|---|---|"]
    for k in sorted(table):
        t = table[k]
        ms = "/".join(f"{t[s]['mix_shift']['power']:.2f}" for s in KEYS)
        mb = "/".join(f"{t[s]['mix_broaden']['power']:.2f}" for s in KEYS)
        fw = max(t[s]["null"]["any"] for s in KEYS)
        fs = max(t[s][sc]["called"] for s in KEYS for sc in t[s] if not sc.startswith("mix_"))
        fdp = np.mean([t[s][sc]["fdp"] for s in KEYS for sc in ("mix_shift", "mix_broaden")])
        bias = np.nanmean([t[s]["mix_shift"]["bias"] for s in KEYS])
        lines.append(f"| {k} | {'yes' if admissible(k) else 'no'} | {ms} | {mb} | {fw:.2f} | {fs:.4f} | {fdp:.3f} | {bias:.1f} |")
    best_e = max((k for k in table if k.startswith("E") and admissible(k)), key=lambda k: power(k, ("mix_shift", "mix_broaden")), default=None)
    best_a = max((f for f in A_GRID if admissible(f"Ac|{f}") and admissible(f"Am|{f}")),
                 key=lambda f: power(f"Ac|{f}", ("mix_shift",)) + power(f"Am|{f}", ("mix_broaden",)), default=None)
    c2 = [k for k in table if k.startswith("C2|") and admissible(k)]
    best_c2 = min(c2, key=lambda k: (round(abs(np.nanmean([table[k][s]["mix_shift"]["bias"] for s in KEYS])), 0),
                                      -power(k, ("mix_shift",))), default=None)
    lines += ["", f"Selected: E = {best_e}; A support fraction = {best_a}; C2 = {best_c2}"]
    (out / "TUNING.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-1:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
