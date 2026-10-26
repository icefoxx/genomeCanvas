"""Validation of the v2.1 per-region shape tests, end to end through ``per_region_tests``.

Fresh seeds (``SEED_BASE``), never used for tuning. Source hash of the profileseg package is recorded
before running; the criteria below were fixed before the first run.

Sources × designs: parametric / real templates / semi-real (real replicate noise) × 3v3 / 2v2.
Scenarios: null; amplitude-only (whole set), site-specific amplitude (parametric/templates) or
global scale (semi-real); 20 % mixtures of shift (100 bp and 50 bp) and broadening (×1.5).

Criteria, per source × design cell (20 runs):
  V1  shape_kernel complete-null FWER ≤ 0.10
  V2  shape_kernel and shift: share of regions called ≤ 0.005 in every no-shape scenario
  V3  shape_kernel and shift: mean FDP ≤ 0.10 in the mixtures
  V4  shift: median estimate over truly shifted regions within ±20 % of the true shift (100 bp mixture)
  Power is reported, not judged.

    python benchmarks/profileseg_region_validation.py --out DIR --templates T.npz --semireal CACHE.npz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

import genomecanvas.profileseg as PS
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, Profiles
from genomecanvas.profileseg.regional import bh, per_region_tests
from genomecanvas.profileseg.simulate import SimSettings, simulate
from semireal import SemiReal

SOURCES = ("param", "templates", "semireal")
DESIGNS = (3, 2)
SCENARIOS = ("null", "amplitude", "other_amplitude", "mix_shift", "mix_shift50", "mix_broaden")
CRIT = dict(V1=0.10, V2=0.005, V3=0.10, V4=0.20)
FRACTION = 0.2
# Round 1 (v2.1): 300 000. Round 2 (v2.2: shift_score calls, Helmert semi-real noise): 700 000.
SEED_BASE = 700_000
_SEMI = {}


def source_hash() -> str:
    h = hashlib.sha256()
    for f in sorted(Path(PS.__file__).parent.glob("*.py")):
        h.update(f.read_bytes())
    return h.hexdigest()[:16]


def _mix(a: Profiles, b: Profiles, k: int, n: int) -> Profiles:
    return Profiles(np.concatenate([a.values[:, :k], b.values[:, : n - k]], axis=1), a.edges, a.samples,
                    a.regions[:k] + b.regions[: n - k], True)


def make(source, n, scenario, seed, templates_path, cache) -> tuple[Profiles, np.ndarray, float]:
    shift = 50.0 if scenario == "mix_shift50" else 100.0
    if source == "semireal":
        if cache not in _SEMI:
            _SEMI[cache] = SemiReal(Path(cache))
        change = {"null": "null", "amplitude": "amplitude", "other_amplitude": "global_scale",
                  "mix_shift": "shift", "mix_shift50": "shift", "mix_broaden": "broaden"}[scenario]
        frac = 1.0 if scenario == "amplitude" else FRACTION
        prof, changed = _SEMI[cache].experiment(change, 300, frac, seed, n_per_group=n, shift_bp=shift)
        return prof, changed, shift
    templates = np.load(templates_path)["templates"] if source == "templates" else None
    st = SimSettings(n_regions=300, n_ref=n, n_alt=n, shift_bp=shift)
    if not scenario.startswith("mix_"):
        name = {"other_amplitude": "site_amplitude"}.get(scenario, scenario)
        prof, _ = simulate(name, st, seed=seed, templates=templates)
        return prof, np.zeros(len(prof.regions), bool), shift
    change = "broaden" if scenario == "mix_broaden" else "shift"
    a, _ = simulate(change, st, seed=seed, templates=templates)
    b, _ = simulate("null", st, seed=seed + 77_777, templates=templates)
    k = int(FRACTION * 300)
    return _mix(a, b, k, 300), np.r_[np.ones(k, bool), np.zeros(300 - k, bool)], shift


def run(job):
    source, n, scenario, rep, templates_path, cache = job
    seed = SEED_BASE + 10_000 * SOURCES.index(source) + 1000 * SCENARIOS.index(scenario) + 100 * n + rep
    prof, changed, shift = make(source, n, scenario, seed, templates_path, cache)
    a = analyse(prof, Config(reference="WT", alternative="KO"), bootstrap=0)
    reg = per_region_tests(a)
    out = {}
    for m in ("shape_kernel", "shift", "centroid", "mad_width", "log2_scale"):
        if m not in reg.tests:
            continue
        t = reg.tests[m]
        idx = t["region_index"]
        truth = changed[idx]
        called = np.isin(t["call"], ("Higher", "Lower", "Changed"))
        nc = int(called.sum())
        est = np.asarray(t["estimate"], float)
        out[m] = dict(called=float(called.mean()), any=float(called.any()),
                      fdp=float((called & ~truth).sum() / nc) if nc else 0.0,
                      power=float((called & truth).sum() / truth.sum()) if truth.any() else float("nan"),
                      est=float(np.nanmedian(est[truth])) if truth.any() and np.isfinite(est).any() else float("nan"))
    return source, n, scenario, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--templates", required=True)
    ap.add_argument("--semireal", required=True)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    sha = source_hash()
    jobs = [(s, n, sc, r, args.templates, args.semireal) for s in SOURCES for n in DESIGNS for sc in SCENARIOS
            for r in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for s, n, sc, res in pool.map(run, jobs):
            for m, v in res.items():
                agg.setdefault((s, n), {}).setdefault(sc, {}).setdefault(m, []).append(v)
    mean = lambda rows, k: float(np.nanmean([r[k] for r in rows]))
    lines = [f"# Per-region shape validation", "", f"Source hash `{sha}`; {args.reps} runs per cell; seeds {SEED_BASE}+.",
             "", "| source | design | V1 null FWER (kernel) | V2 max false share (kernel / shift) | V3 mean FDP (kernel / shift) | "
             "V4 shift estimate (100 bp) | pass | power kernel: shift100 / shift50 / broaden | power shift: 100 / 50 |",
             "|---|---|---|---|---|---|---|---|---|"]
    verdict = {}
    for (s, n), d in sorted(agg.items()):
        v1 = mean(d["null"]["shape_kernel"], "any")
        no_shape = ("null", "amplitude", "other_amplitude")
        v2k = max(mean(d[sc]["shape_kernel"], "called") for sc in no_shape)
        v2s = max(mean(d[sc]["shift"], "called") for sc in no_shape)
        mixes = ("mix_shift", "mix_shift50", "mix_broaden")
        v3k = max(mean(d[sc]["shape_kernel"], "fdp") for sc in mixes)
        v3s = max(mean(d[sc]["shift"], "fdp") for sc in mixes)
        v4 = mean(d["mix_shift"]["shift"], "est")
        ok = dict(V1=v1 <= CRIT["V1"], V2=max(v2k, v2s) <= CRIT["V2"], V3=max(v3k, v3s) <= CRIT["V3"],
                  V4=abs(v4 - 100) <= CRIT["V4"] * 100)
        verdict[f"{s}:{n}v{n}"] = dict(V1=v1, V2=[v2k, v2s], V3=[v3k, v3s], V4=v4, passed=ok)
        pk = "/".join(f"{mean(d[sc]['shape_kernel'], 'power'):.2f}" for sc in mixes)
        ps = "/".join(f"{mean(d[sc]['shift'], 'power'):.2f}" for sc in ("mix_shift", "mix_shift50"))
        fails = [k for k, v in ok.items() if not v]
        lines.append(f"| {s} | {n}v{n} | {v1:.2f} | {v2k:.4f} / {v2s:.4f} | {v3k:.3f} / {v3s:.3f} | {v4:.1f} | "
                     f"{'all' if not fails else 'FAIL ' + ','.join(fails)} | {pk} | {ps} |")
    (out / "region_validation.json").write_text(json.dumps(dict(sha=sha, criteria=CRIT, cells=verdict), indent=1))
    (out / "REGION_VALIDATION.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
