"""One replicate with a sample-wide artefact: effect on per-region shape calls, and detection.

A sample-wide artefact (one replicate shifted by 30 bp in every region, or broadened ×1.25 as a
different fragment length would do) is injected into the first alternative-group sample.
Scenarios: null, null + artefact, 20 % shifted (100 bp), 20 % shifted + artefact.

Reported: false-call share / FWER / FDP / power of `shape_kernel` and `shift`, and how often the
replicate-influence diagnostics flag the run (sample-wide offset flag; leave-one-out flag).

    python benchmarks/profileseg_replicate_artefact.py --out DIR --templates T.npz --semireal CACHE.npz
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

from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, Profiles
from genomecanvas.profileseg.regional import per_region_tests
from profileseg_region_validation import make

SOURCES = ("param", "templates", "semireal")
DESIGNS = (3, 2)
SCENARIOS = ("null", "null+shift_artefact", "null+width_artefact", "mix", "mix+shift_artefact")


def inject(prof: Profiles, sample: int, kind: str) -> Profiles:
    v = prof.values.copy()
    c = 0.5 * (prof.edges[:-1] + prof.edges[1:])
    for r in range(v.shape[1]):
        x = np.nan_to_num(v[sample, r])
        if kind == "shift":
            v[sample, r] = np.interp(c - 30.0, c, x, left=x[0], right=x[-1])
        else:
            bg = np.r_[x[:25], x[-25:]].mean()
            y = bg + np.interp(c / 1.25, c, x - bg)
            v[sample, r] = y * max(x.sum(), 1e-9) / max(y.sum(), 1e-9)
    return Profiles(v, prof.edges, prof.samples, prof.regions, prof.oriented)


def run(job):
    source, n, scenario, rep, tpl, cache = job
    seed = 900_000 + 10_000 * SOURCES.index(source) + 1000 * SCENARIOS.index(scenario) + 100 * n + rep
    base = "mix_shift" if scenario.startswith("mix") else "null"
    prof, changed, _ = make(source, n, base, seed, tpl, cache)
    if "artefact" in scenario:
        prof = inject(prof, n, "shift" if "shift_artefact" in scenario else "width")  # first alternative sample
    a = analyse(prof, Config(reference="WT", alternative="KO"), bootstrap=0)
    reg = per_region_tests(a)
    out = {}
    for m in ("shape_kernel", "shift"):
        t = reg.tests[m]
        truth = changed[t["region_index"]]
        called = np.isin(t["call"], ("Higher", "Lower", "Changed"))
        nc = int(called.sum())
        out[m] = dict(called=float(called.mean()), any=float(called.any()),
                      fdp=float((called & ~truth).sum() / nc) if nc else 0.0,
                      power=float((called & truth).sum() / truth.sum()) if truth.any() else float("nan"))
    flags = reg.influence.get("flags", [])
    out["flag_offset"] = float(any("sample-wide" in f for f in flags))
    out["flag_offset_right_sample"] = float(any("sample-wide" in f and "KO_rep1" in f for f in flags) and not any("sample-wide" in f and "KO_rep1" not in f for f in flags))
    out["flag_loo"] = float(any("depend on that replicate" in f for f in flags))
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
            agg.setdefault(f"{s}:{n}v{n}", {}).setdefault(sc, []).append(res)
    lines = ["# One replicate with a sample-wide artefact", "",
             "| cell | scenario | kernel: called / FWER / FDP / power | shift: called / FWER / FDP / power | "
             "offset flag (right sample) | leave-one-out flag |", "|---|---|---|---|---|---|"]
    table = {}
    for cell in sorted(agg):
        for sc in SCENARIOS:
            rows = agg[cell][sc]
            m = lambda k1, k2: float(np.nanmean([r[k1][k2] if isinstance(r[k1], dict) else r[k1] for r in rows]))
            f = lambda k: float(np.mean([r[k] for r in rows]))
            table.setdefault(cell, {})[sc] = {k: (m(k, "called"), m(k, "any"), m(k, "fdp"), m(k, "power"))
                                              for k in ("shape_kernel", "shift")}
            table[cell][sc].update(flag_offset=f("flag_offset"), flag_right=f("flag_offset_right_sample"),
                                   flag_loo=f("flag_loo"))
            t = table[cell][sc]
            fmt = lambda x: f"{x[0]:.3f} / {x[1]:.2f} / {x[2]:.3f} / {x[3]:.2f}"
            lines.append(f"| {cell} | {sc} | {fmt(t['shape_kernel'])} | {fmt(t['shift'])} | "
                         f"{t['flag_offset']:.2f} ({t['flag_right']:.2f}) | {t['flag_loo']:.2f} |")
    (out / "replicate_artefact.json").write_text(json.dumps(table, indent=1))
    (out / "REPLICATE_ARTEFACT.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
