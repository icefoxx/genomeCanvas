"""Calibration boundaries of the per-region shape tests (pooled permutation null), at nominal FDR 5 %.

Design and criteria are fixed here, before any result is seen (seeds 5 000 000+, not used before).

Baseline: parametric shapes, 600 regions, 20 % changed (half translated by 100 bp, half broadened
×1.5), 2 vs 2 and 3 vs 3, 100 runs per cell (raised from 30 after a one-run smoke test, for precision only; criteria unchanged). One factor is varied at a time:

  width        peak sd 60 / 150 / 400 bp
  depth        reads per bin at the apex 5 / 15 / 40 / 160
  noise        low (replicate sd 0.05, AR sd 0.04, dispersion 0.05) / default / high (0.3, 0.15, 0.4)
  proportion   0 % (complete null) / 5 % / 20 % / 50 % of regions changed
  correlation  every region duplicated with its anchor offset by 200 bp (overlapping windows)
  mixed        one run mixes three region classes (narrow 60 bp, broad 400 bp, low depth 5):
               calibration reported per class
  batch        one replicate per group carries a sample-wide 15 bp translation artefact (balanced)
  semireal     real replicate noise (GSE112560 Helmert residuals), 300 regions
  mixed_large  (protocol v2, added after the first run) the mixed design with 1000 regions per class,
               where the covariate-stratified null (v3 default) has enough regions to stratify
  mixed_n      (protocol v3, review round 3) the mixed design at 600 / 1000 / 1500 / 2000 / 3000 / 6000
               regions in total, to locate the region count at which per-class calibration holds

Outcomes per test (shape_kernel, shift, shape_any): FDP (mean and 95 % Monte Carlo CI over runs),
false-call share among unchanged regions, share of unchanged regions with p < 0.01 (per-class tail
calibration), complete-null FWER (binomial CI), power.

Criteria (nominal FDR 0.05):
  EXCEEDS   the 95 % CI of the mean FDP lies above 0.05 (or, at 0 %, the FWER CI lies above 0.05)
  WATCH     point estimate above 0.075 without the CI excluding 0.05
  ok        otherwise
"ok" means no statistically clear exceedance in that cell; it does not establish FDR ≤ 0.05.

Options: --null strata|none (v3 / v2.4 pooled null), --fdr-strata none|spread|depth, --seed-base,
--factors (subset). Every output records the git commit, uncommitted changes under src/, options,
seeds and run times; the Markdown table is generated from the JSON written in the same run.

    python benchmarks/profileseg_calibration_boundary.py --out DIR --semireal CACHE.npz [--reps 100]
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import platform
import subprocess
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from scipy import stats

from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, Profiles
from genomecanvas.profileseg.regional import per_region_tests
from genomecanvas.profileseg.simulate import SimSettings, simulate

TESTS = ("shape_kernel", "shift", "shape_any")
NOISE = {"low": dict(replicate_sd=0.05, ar_sd=0.04, dispersion=0.05), "default": {},
         "high": dict(replicate_sd=0.3, ar_sd=0.15, dispersion=0.4)}
CELLS = ([("width", w) for w in (60, 150, 400)] + [("depth", d) for d in (5, 15, 40, 160)]
         + [("noise", k) for k in NOISE] + [("proportion", p) for p in (0.0, 0.05, 0.2, 0.5)]
         + [("correlation", "overlap"), ("mixed", "3 classes"), ("batch", "balanced 15 bp"), ("semireal", "real noise"),
            ("mixed_large", "3 classes × 1000")]
         + [("mixed_n", t) for t in (600, 1000, 1500, 2000, 3000, 6000)])
SEED = 5_000_000


def _concat(parts: list[Profiles]) -> Profiles:
    values = np.concatenate([p.values for p in parts], axis=1)
    regions = [dict(chrom="c", start=10_000 * (i + 1), end=10_000 * (i + 1) + 1, name=f"r{i}", strand="+")
               for i in range(values.shape[1])]
    frags = (np.concatenate([p.fragments for p in parts], axis=1)
             if all(getattr(p, "fragments", None) is not None for p in parts) else None)
    return Profiles(values, parts[0].edges, parts[0].samples, regions, True, fragments=frags)


def _part(name: str, n: int, k: int, seed: int, **kw) -> Profiles:
    """Simulated part with its count depth attached (counts per region and sample, as BAM input gives)."""
    if not k:
        return None
    prof, truth = simulate(name, SimSettings(n_regions=k, n_ref=n, n_alt=n, **kw), seed=seed)
    lib = np.asarray(truth["library"], float)
    prof.fragments = np.rint(np.nan_to_num(prof.values) * lib[:, None, None]).sum(-1)
    return prof


def build(factor: str, level, n: int, seed: int, semireal: str | None):
    """Profiles, truth (changed), class labels."""
    N, frac, kw = 600, 0.2, {}
    if factor == "width":
        kw = dict(peak_sigma=float(level))
    elif factor == "depth":
        kw = dict(depth=float(level))
    elif factor == "noise":
        kw = NOISE[level]
    elif factor == "proportion":
        frac = float(level)
    if factor == "semireal":
        import sys
        sys.path.insert(0, str(Path(__file__).parent))
        from semireal import SemiReal
        sr = SemiReal(Path(semireal))
        a, ch_a = sr.experiment("shift", 150, 0.4, seed, n_per_group=n)   # 60 shifted of 150
        b, ch_b = sr.experiment("broaden", 150, 0.4, seed + 1, n_per_group=n)
        c, _ = sr.experiment("null", 300, 0.0, seed + 2, n_per_group=n)
        prof = _concat([a, b, c])
        changed = np.r_[ch_a, ch_b, np.zeros(300, bool)]
        return prof, changed, np.zeros(len(changed), int)
    if factor in ("mixed", "mixed_large", "mixed_n"):
        parts, changed, cls = [], [], []
        for j, extra in enumerate((dict(peak_sigma=60.0), dict(peak_sigma=400.0), dict(depth=5.0))):
            k = 1000 if factor == "mixed_large" else int(level) // 3 if factor == "mixed_n" else 200
            m = int(round(frac * k))
            for name, cnt in (("shift", m // 2), ("broaden", m - m // 2), ("null", k - m)):
                if cnt:
                    parts.append(_part(name, n, cnt, seed + 31 * j + zlib.crc32(name.encode()) % 7, **extra))
                    changed += [name != "null"] * cnt
                    cls += [j] * cnt
        return _concat(parts), np.array(changed), np.array(cls)
    m = int(round(frac * N))
    parts, changed = [], []
    for j, (name, cnt) in enumerate((("shift", m // 2), ("broaden", m - m // 2), ("null", N - m))):
        if cnt:
            parts.append(_part(name, n, cnt, seed + 101 * j, **kw))
            changed += [name != "null"] * cnt
    prof = _concat(parts)
    changed = np.array(changed)
    if factor == "correlation":  # duplicate every region with its anchor moved by 200 bp
        sh = int(200 / float(np.median(np.diff(prof.edges))))
        dup = np.concatenate([prof.values[..., sh:], np.repeat(prof.values[..., -1:], sh, axis=-1)], axis=-1)
        prof = _concat([prof, Profiles(dup, prof.edges, prof.samples, prof.regions, True, fragments=prof.fragments)])
        changed = np.r_[changed, changed]
    if factor == "batch":  # sample-wide 15 bp translation in the first replicate of each group
        c = 0.5 * (prof.edges[:-1] + prof.edges[1:])
        v = prof.values.copy()
        for s in (0, n):
            for r in range(v.shape[1]):
                x = np.nan_to_num(v[s, r])
                v[s, r] = np.interp(c - 15.0, c, x, left=x[0], right=x[-1])
        prof = Profiles(v, prof.edges, prof.samples, prof.regions, True, fragments=prof.fragments)
    return prof, changed, np.zeros(len(changed), int)


def run(job):
    factor, level, n, rep, semireal, null, fdr_strata, seed_base = job
    import genomecanvas.profileseg.regional as R
    R.NULL_ADJUST = null
    seed = seed_base + 10_000 * [c[0] for c in CELLS].index(factor) + 1000 * rep + 100 * n + zlib.crc32(str(level).encode()) % 97
    prof, changed, cls = build(factor, level, n, seed, semireal)
    a = analyse(prof, Config(reference="WT", alternative="KO", fdr_strata=fdr_strata), bootstrap=0)
    reg = per_region_tests(a)
    out = {}
    for m in TESTS:
        t = reg.tests[m]
        idx = t["region_index"]
        called = np.isin(t["call"], ("Higher", "Lower", "Changed"))
        tr, cl = changed[idx], cls[idx]
        res = {}
        for label, sel in [("all", np.ones(len(idx), bool))] + [(f"class{j}", cl == j) for j in np.unique(cl) if len(np.unique(cl)) > 1]:
            c, t_ = called[sel], tr[sel]
            nc = int(c.sum())
            pv = np.asarray(t["p"], float)[sel][~t_]
            res[label] = dict(fdp=float((c & ~t_).sum() / nc) if nc else 0.0, any=float((c & ~t_).any()),
                              p01=float(np.mean(pv[np.isfinite(pv)] < 0.01)) if np.isfinite(pv).any() else float("nan"),
                              false_share=float((c & ~t_).sum() / max((~t_).sum(), 1)),
                              power=float((c & t_).sum() / t_.sum()) if t_.any() else float("nan"))
        out[m] = res
    return factor, str(level), n, out


def ci_mean(x):
    x = np.asarray(x, float)
    m = float(np.mean(x))
    h = 1.96 * float(np.std(x, ddof=1)) / np.sqrt(len(x)) if len(x) > 1 else float("nan")
    return m, m - h, m + h


def ci_binom(k, n):
    lo = stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return k / n, float(lo), float(hi)


def verdict(m, lo, hi):
    if lo > 0.05:
        return "EXCEEDS"
    if m > 0.075:
        return "WATCH"
    return "ok"


def provenance(args) -> dict:
    root = Path(__file__).resolve().parents[1]
    git = lambda *a: subprocess.run(["git", *a], cwd=root, capture_output=True, text=True).stdout.strip()  # noqa: E731
    reg = root / "src/genomecanvas/profileseg/regional.py"
    return dict(commit=git("rev-parse", "HEAD"), branch=git("rev-parse", "--abbrev-ref", "HEAD"),
                uncommitted_src=[ln for ln in git("status", "--porcelain", "--", "src").splitlines()],
                regional_py_sha256=hashlib.sha256(reg.read_bytes()).hexdigest()[:16],
                src_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--", "src").encode()).hexdigest()[:16],
                python=platform.python_version(), numpy=np.__version__,
                options=dict(reps=args.reps, null=args.null, fdr_strata=args.fdr_strata, seed_base=args.seed_base,
                             factors=args.factors or "all"))


def summarise(agg: dict) -> dict:
    """Cell summaries (estimates and CIs) from per-run outcomes."""
    out = {}
    for (f, lv, n), rows in agg.items():
        for label in rows[0][TESTS[0]]:
            cell = {}
            for t in TESTS:
                r_ = [r[t][label] for r in rows]
                with np.errstate(all="ignore"):
                    pw = float(np.nanmean([x["power"] for x in r_])) if any(np.isfinite(x["power"]) for x in r_) else float("nan")
                d = dict(runs=len(r_), power=pw,
                         p01_unchanged=float(np.nanmean([x["p01"] for x in r_])),
                         false_share=float(np.mean([x["false_share"] for x in r_])))
                if f == "proportion" and float(lv) == 0.0:
                    k = int(sum(x["any"] for x in r_))
                    d["fwer"], d["low"], d["high"] = ci_binom(k, len(r_))
                    d["verdict"] = verdict(d["fwer"], d["low"], d["high"])
                else:
                    d["fdp"], d["low"], d["high"] = ci_mean([x["fdp"] for x in r_])
                    d["verdict"] = verdict(d["fdp"], d["low"], d["high"])
                cell[t] = d
            out[f"{f}|{lv}|{n}v{n}|{label}"] = cell
    return out


def to_markdown(doc: dict) -> str:
    pv = doc["provenance"]
    lines = ["# Calibration boundaries of the per-region shape tests (nominal FDR 0.05)", "",
             f"Generated from `calibration_boundary.json` of the same run. Commit `{pv['commit'][:10]}` "
             f"({pv['branch']}); uncommitted changes under src/: {len(pv['uncommitted_src'])}; regional.py "
             f"sha256 `{pv['regional_py_sha256']}`; src diff sha256 `{pv['src_diff_sha256']}`; options {json.dumps(pv['options'])}; "
             f"run {doc['started']} → {doc['finished']}.", "",
             "FDP: mean [95 % Monte Carlo CI] / power; at 0 % changed the FWER [exact binomial CI]. "
             "Verdict: EXCEEDS = CI above 0.05; WATCH = estimate above 0.075; ok = no clear exceedance "
             "(not proof of FDR ≤ 0.05). p01 = share of unchanged regions with p < 0.01 (nominal 0.01).", "",
             "| factor | level | design | class | " + " | ".join(f"{t}: FDP / power / p01" for t in TESTS) + " | verdict |",
             "|---|---|---|---|" + "---|" * (len(TESTS) + 1)]
    for key, cell in doc["cells"].items():
        f, lv, des, label = key.split("|")
        parts, bad = [], []
        for t in TESTS:
            d = cell[t]
            if "fwer" in d:
                parts.append(f"FWER {d['fwer']:.3f} [{d['low']:.3f}, {d['high']:.3f}] / – / {d['p01_unchanged']:.4f}")
            else:
                parts.append(f"{d['fdp']:.3f} [{max(d['low'], 0):.3f}, {d['high']:.3f}] / {d['power']:.2f} / {d['p01_unchanged']:.4f}")
            if d["verdict"] != "ok":
                bad.append(f"{t}:{d['verdict']}")
        lines.append(f"| {f} | {lv} | {des} | {label} | " + " | ".join(parts) + f" | {'; '.join(bad) or 'ok'} |")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--semireal", required=True)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--null", choices=("strata", "none"), default="strata")
    ap.add_argument("--fdr-strata", choices=("none", "spread", "depth"), default="none")
    ap.add_argument("--seed-base", type=int, default=SEED)
    ap.add_argument("--factors", default="", help="comma-separated subset of factors")
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    cells = [c for c in CELLS if not args.factors or c[0] in args.factors.split(",")]
    started = datetime.datetime.now().isoformat(timespec="seconds")
    prov = provenance(args)  # read before the run: the code that runs, not the code at the end
    jobs = [(f, lv, n, r, args.semireal, args.null, args.fdr_strata, args.seed_base)
            for f, lv in cells for n in (2, 3) for r in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(args.jobs) as ex:
        for f, lv, n, res in ex.map(run, jobs):
            agg.setdefault((f, lv, n), []).append(res)
    doc = dict(provenance=prov, started=started,
               finished=datetime.datetime.now().isoformat(timespec="seconds"), cells=summarise(agg))
    (out / "calibration_boundary.json").write_text(json.dumps(doc, indent=1))
    md = to_markdown(doc)
    (out / "CALIBRATION_BOUNDARY.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
