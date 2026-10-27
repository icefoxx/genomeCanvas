"""Shift size × peak shape: power of the shape tests and accuracy of the translation estimate.

20 % of 600 regions are translated by d bp (d = 10 … 200), the rest unchanged; 2 vs 2 and 3 vs 3,
30 runs per cell. Shapes: Gaussian peaks (sd 60 / 150 / 400 bp), bimodal templates (two peaks;
profileseg_validation.py) and real templates (GSE112560, if given). Reported per test (`shift`,
`shape_kernel`, `shape_any`): power and FDP at BH 5 %; for the translation estimate, the median of
estimate / true shift over all translated regions and over called ones.

The `shift` test projects each profile on the derivative of the pooled template: locally most
powerful for small translations, and expected to lose power once d is comparable to the peak width
or, for periodic profiles, to half the spacing. This script measures where.

v3 translation estimate (cross-correlation of condition means, with reliability flags): reported per
cell are the share of called regions flagged reliable, the direction accuracy among them, and the
median estimate ÷ true over reliable called regions, next to the template-fit estimate. A non-translation
control ("side_gain": a one-sided added peak, half the apex, d bp downstream of the summit) checks how
often a shape change that is not a translation is flagged as a reliable translation.

Criteria fixed before the v3 run (seeds 9 400 000+): direction accuracy among reliable calls ≥ 0.95;
median reliable estimate ÷ true within [0.85, 1.15] for d ≥ 50 bp; side_gain: ≤ 0.2 of called regions
flagged reliable. Round 5 (translation_r2 ≥ 0.8, chosen on development seeds 7 700 000+): the same three
criteria on new seeds (--seed-base 9700000), and the share of called translations kept as consistent reported.

    python benchmarks/profileseg_shift_boundary.py --out DIR [--templates T.npz] [--reps 30]
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.regional import per_region_tests
from genomecanvas.profileseg.simulate import SimSettings, simulate

SHIFTS = (10, 20, 35, 50, 75, 100, 150, 200)
TESTS = ("shift", "shape_kernel", "shape_any")
SEED = 9_400_000
MODES = ("translate", "side_gain")


def shapes(templates_path: str | None) -> dict:
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from profileseg_validation import bimodal_templates
    base = SimSettings()
    out = {f"gauss{s}": dict(peak_sigma=float(s)) for s in (60, 150, 400)}
    out["bimodal"] = dict(templates=bimodal_templates(1000, base.flank, base.binsize, 11))
    if templates_path:
        out["real_GSE112560"] = dict(templates=np.load(templates_path)["templates"])
    return out


def run(job):
    shape, kw, d, n, rep, mode = job[:6]
    base = job[6] if len(job) > 6 else SEED
    kw = dict(kw)
    templates = kw.pop("templates", None)
    seed = (base + 100_000 * n + 1000 * SHIFTS.index(d) + rep + 10_000 * (zlib.crc32(shape.encode()) % 50)
            + 5_000_000 * MODES.index(mode))
    st = SimSettings(n_regions=600, n_ref=n, n_alt=n, shift_bp=float(d), sparse_fraction=0.2, **kw)
    if mode == "translate":
        prof, truth = simulate("sparse_shift", st, seed=seed, templates=templates)
        changed = np.asarray(truth["regions"]["shift"]) != 0
    else:  # one-sided gain: a peak of half the apex added d bp downstream of the summit, alternative only
        prof, truth = simulate("null", st, seed=seed, templates=templates)
        rng = np.random.default_rng(seed)
        changed = np.zeros(prof.values.shape[1], bool)
        changed[rng.choice(len(changed), int(0.2 * len(changed)), replace=False)] = True
        c = 0.5 * (prof.edges[:-1] + prof.edges[1:])
        alt = np.flatnonzero(prof.conditions == "KO")
        v = prof.values.copy()
        mean = np.nan_to_num(v).mean(0)
        for r in np.flatnonzero(changed):
            bg = np.median(np.r_[mean[r, :10], mean[r, -10:]])
            k = int(np.argmax(mean[r]))
            bump = 0.5 * max(mean[r, k] - bg, 0) * np.exp(-0.5 * ((c - c[k] - d) / 40.0) ** 2)
            v[alt, r] += bump
        prof.values = v
    reg = per_region_tests(analyse(prof, Config(reference="WT", alternative="KO"), bootstrap=0))
    out = {}
    for m in TESTS:
        t = reg.tests[m]
        idx = t["region_index"]
        c = np.isin(t["call"], ("Higher", "Lower", "Changed"))
        tr = changed[idx]
        nc = int(c.sum())
        out[m] = dict(power=float((c & tr).sum() / max(tr.sum(), 1)), fdp=float((c & ~tr).sum() / nc) if nc else 0.0)
        if m == "shift":
            est = np.asarray(t["estimate"], float)  # NaN where not reliable
            fit = np.asarray(t["shift_bp_fit"], float)
            rel = np.asarray(t["translation_consistent"], bool)
            ct = c & tr
            nan = float("nan")
            out[m]["reliable_share_called"] = float(rel[ct].mean()) if ct.any() else nan
            rc = ct & rel
            out[m]["direction_reliable"] = float(np.mean(est[rc] > 0)) if rc.any() else nan
            out[m]["label_direction"] = float(np.mean(t["call"][rc] == "Higher")) if rc.any() else nan
            out[m]["ratio_reliable"] = float(np.median(est[rc] / d)) if rc.any() else nan
            out[m]["fit_ratio_called"] = float(np.median(fit[ct] / d)) if ct.any() else nan
            out[m]["fit_direction_called"] = float(np.mean(fit[ct] > 0)) if ct.any() else nan
    return shape, d, n, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--templates")
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--seed-base", type=int, default=SEED)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    git = lambda *a: subprocess.run(["git", *a], cwd=root, capture_output=True, text=True).stdout.strip()  # noqa: E731
    started = datetime.datetime.now().isoformat(timespec="seconds")
    prov = dict(  # read before the run
        commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                regional_py_sha256=hashlib.sha256((root / "src/genomecanvas/profileseg/regional.py").read_bytes()).hexdigest()[:16],
                src_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--", "src").encode()).hexdigest()[:16],
                options=dict(reps=args.reps, seed_base=args.seed_base, templates=bool(args.templates)))
    sh = shapes(args.templates)
    jobs = [(s, kw, d, n, r, mode, args.seed_base) for mode in MODES for s, kw in sh.items() for d in SHIFTS for n in (2, 3)
            for r in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(args.jobs) as ex:
        for job, (s, d, n, res) in zip(jobs, ex.map(run, jobs), strict=True):
            agg.setdefault((job[5], s, d, n), []).append(res)
    cells = {}
    for (mode, s, d, n), rows in agg.items():
        c = {}
        for m in TESTS:
            with np.errstate(all="ignore"):
                c[m] = {k: (float(np.nanmean(v)) if np.isfinite(v := np.array([r[m][k] for r in rows], float)).any() else float("nan"))
                        for k in rows[0][m]}
        cells[f"{mode}|{s}|{d}|{n}v{n}"] = c
    doc = dict(provenance=prov, started=started, finished=datetime.datetime.now().isoformat(timespec="seconds"), cells=cells)
    (out / "shift_boundary.json").write_text(json.dumps(doc, indent=1))
    lines = ["# Shift size × peak shape (20 % of 600 regions changed; BH 5 %)", "",
             f"Commit `{prov['commit'][:10]}`; uncommitted src files {prov['uncommitted_src']}; regional.py `{prov['regional_py_sha256']}`; "
             f"src diff `{prov['src_diff_sha256']}`; {args.reps} runs per cell; seeds {args.seed_base}+; {started} → {doc['finished']}.", "",
             "Power and FDP per test; for `shift` among called changed regions: share flagged reliable, direction accuracy "
             "of the reported estimate among reliable ones (the Higher/Lower label is shown after the slash), median reported "
             "estimate ÷ true among reliable ones, and the template fit (median ÷ true, sign accuracy) for comparison.", "",
             "| mode | shape | design | d (bp) | power shift / kernel / any | FDP shift / kernel / any | reliable share | direction (estimate / label) | reliable estimate ÷ true | template fit ÷ true (sign ok) |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for mode in MODES:
        for s in sh:
            for n in (2, 3):
                for d in SHIFTS:
                    c = cells[f"{mode}|{s}|{d}|{n}v{n}"]
                    x = c["shift"]
                    lines.append(f"| {mode} | {s} | {n}v{n} | {d} | {x['power']:.2f} / {c['shape_kernel']['power']:.2f} / {c['shape_any']['power']:.2f} | "
                                 f"{x['fdp']:.3f} / {c['shape_kernel']['fdp']:.3f} / {c['shape_any']['fdp']:.3f} | {x['reliable_share_called']:.2f} | "
                                 f"{x['direction_reliable']:.3f} / {x['label_direction']:.3f} | {x['ratio_reliable']:.2f} | "
                                 f"{x['fit_ratio_called']:.2f} ({x['fit_direction_called']:.2f}) |")
    (out / "SHIFT_BOUNDARY.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
