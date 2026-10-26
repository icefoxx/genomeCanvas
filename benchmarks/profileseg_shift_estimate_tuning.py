"""Tune the `shift` effect-size estimator (template fit): template smoothing and condition-blind
template refinement. Tuning seeds (2 100 000+); calls are unaffected (they come from the score test).

Reports, over truly shifted regions, the median estimate relative to the true shift and the median
absolute error, for 20 % mixtures with 100 bp and 50 bp shifts.

    python benchmarks/profileseg_shift_estimate_tuning.py --templates T.npz --semireal CACHE.npz
"""

from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.regional import template_fit
from profileseg_region_validation import make

SOURCES = ("param", "templates", "semireal")
VARIANTS = [(sm, rf) for sm in (1, 3, 5) for rf in (0, 1, 2)] + [(7, 2), (9, 2), (7, 3), (11, 2)]


def run(job):
    source, n, scenario, rep, tpl, cache = job
    seed = 2_100_000 + 10_000 * SOURCES.index(source) + 1000 * (scenario == "mix_shift50") + 100 * n + rep
    prof, changed, true_shift = make(source, n, scenario, seed, tpl, cache)
    cfg = Config(reference="WT", alternative="KO")
    a = analyse(prof, cfg, bootstrap=0)
    comp = np.flatnonzero(a.complete)
    sidx = comp[a.core_mass[:, comp].mean(0) > 0]
    y = a.core_signal[:, sidx]
    core = M.core_mask(prof.edges, cfg.core_flank)
    cov = np.nan_to_num(a.values[:, sidx][..., core])
    truth = changed[sidx]
    out = {}
    for sm, rf in VARIANTS:
        f = template_fit(y, cov, a.core_edges, smooth=sm, refine=rf)
        est = f["shift"][n:].mean(0) - f["shift"][:n].mean(0)
        e = est[truth]
        out[(sm, rf)] = (float(np.median(e) / true_shift), float(np.median(np.abs(e - true_shift)) / true_shift))
    return source, n, scenario, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--templates", required=True)
    ap.add_argument("--semireal", required=True)
    ap.add_argument("--reps", type=int, default=6)
    args = ap.parse_args(argv)
    jobs = [(s, n, sc, r, args.templates, args.semireal) for s in SOURCES for n in (2, 3)
            for sc in ("mix_shift", "mix_shift50") for r in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(6) as ex:
        for s, n, sc, out in ex.map(run, jobs):
            for v, x in out.items():
                agg.setdefault(v, {}).setdefault((s, n, sc), []).append(x)
    cells = sorted(next(iter(agg.values())))
    print("| smooth, refine | " + " | ".join(f"{s} {n}v{n} {'100' if sc == 'mix_shift' else '50'} bp" for s, n, sc in cells) + " | worst ratio | mean MAE |")
    print("|---|" + "---|" * (len(cells) + 2))
    for v in VARIANTS:
        r = [np.mean([x[0] for x in agg[v][c]]) for c in cells]
        mae = np.mean([np.mean([x[1] for x in agg[v][c]]) for c in cells])
        worst = min(r, key=lambda z: -abs(z - 1))
        print(f"| {v} | " + " | ".join(f"{x:.2f}" for x in r) + f" | {worst:.2f} | {mae:.2f} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
