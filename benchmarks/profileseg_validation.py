"""Independent validation of ProfileSeg on designs not used during development.

The method is frozen before this runs: the SHA-256 of the ProfileSeg sources is written into the
report, and nothing in the package may be tuned on these results. Acceptance criteria are fixed
here, before any result is seen:

Aggregate (global) inference, per design
  A1  median false-call rate over metrics whose truth is ≈ 0       ≤ 0.07
  A2  maximum false-call rate of any such metric                     ≤ 0.20
  A3  minimum 95 % CI coverage                                       ≥ 0.85
  A4  family-wise false-call rate of the report (primary + Holm)     ≤ 0.10 in every scenario
Per-region inference, per design
  R1  FDR (mean FDP) where true changes exist                        ≤ 0.075 for every metric
  R2  complete-null FWER (= FDR) per metric                          ≤ 0.10
  R3  power of log2_signal_relative in site_amplitude                 reported (no threshold)

Designs
  V1  different noise: dispersion 0.4, replicate sd 0.3, library sd 0.25, AR ρ 0.85 / sd 0.15
  V2a bimodal (nucleosome-like) shapes      V2b right-skewed shapes
  V3a unbalanced 3 vs 2                      V3b 4 vs 4 with 10 bp bins over ±3 kb
  V4  real templates from an independent assay (2-cell IP, GSM5746912/3), if given

Usage:  python benchmarks/profileseg_validation.py --out validation [--reps 40] [--templates-2cell t.npz]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from genomecanvas.profileseg import __file__ as package_init
from genomecanvas.profileseg.benchmark import familywise, run_benchmark, run_regional_benchmark
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.simulate import SimSettings

CRITERIA = dict(A1=0.07, A2=0.20, A3=0.85, A4=0.10, R1=0.075, R2=0.10)


def source_hash() -> str:
    digest = hashlib.sha256()
    for path in sorted(Path(package_init).parent.glob("*.py")):
        digest.update(path.name.encode() + path.read_bytes())
    return digest.hexdigest()


def bimodal_templates(n: int, flank: int, binsize: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    c = np.arange(-flank, flank, binsize) + binsize / 2
    gap = rng.normal(150, 20, n)[:, None]
    left = rng.uniform(0.4, 1.0, n)[:, None]
    t = left * np.exp(-0.5 * ((c + gap) / 60) ** 2) + np.exp(-0.5 * ((c - gap) / 60) ** 2)
    return t / t.max(axis=1, keepdims=True) + 0.3


def skewed_templates(n: int, flank: int, binsize: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    c = np.arange(-flank, flank, binsize) + binsize / 2
    k = rng.uniform(1.5, 3.0, n)[:, None]
    scale = rng.uniform(80, 160, n)[:, None]
    x = np.clip(c[None, :] + 150, 0, None) / scale
    t = np.where(x > 0, x ** (k - 1) * np.exp(-x), 0.0)  # gamma-shaped: steep left, long right tail
    return t / t.max(axis=1, keepdims=True) + 0.3


def designs(two_cell: Path | None) -> dict:
    base = SimSettings()
    out = {
        "V1_noise": dict(settings=SimSettings(dispersion=0.4, replicate_sd=0.3, library_sd=0.25, ar_rho=0.85, ar_sd=0.15)),
        "V2a_bimodal": dict(settings=base, templates=bimodal_templates(1000, base.flank, base.binsize, 11)),
        "V2b_skewed": dict(settings=base, templates=skewed_templates(1000, base.flank, base.binsize, 12)),
        "V3a_3v2": dict(settings=SimSettings(n_alt=2)),
        "V3b_4v4_10bp_3kb": dict(settings=SimSettings(n_ref=4, n_alt=4, flank=3000, binsize=10)),
    }
    if two_cell:
        out["V4_2cell_templates"] = dict(settings=base, templates=np.load(two_cell)["templates"])
    return out


def evaluate(agg: dict, reg: dict) -> dict:
    core = [r for r in agg["summary"] if not r["metric"].startswith("naive")]
    null = [r["welch_flag_rate"] for r in core if r["expected"] == 0]
    cov = [r["ci_coverage"] for r in core if np.isfinite(r["ci_coverage"])]
    fw = familywise(agg["records"])
    with_truth = [r for r in reg["summary"] if r["true_changes"] >= 1]
    complete_null = [r for r in reg["summary"] if r["true_changes"] < 1]
    site = next((r for r in reg["summary"] if r["scenario"] == "site_amplitude" and r["metric"] == "log2_signal_relative"), None)
    values = dict(
        A1=float(np.median(null)), A2=float(np.max(null)), A3=float(np.min(cov)),
        A4=float(max(f["familywise_false_call_rate"] for f in fw)),
        R1=float(max(r["fdr"] for r in with_truth)) if with_truth else float("nan"),
        R2=float(max(r["fdr"] for r in complete_null)) if complete_null else float("nan"),
        R3=float(site["power"]) if site else float("nan"),
    )
    passed = dict(
        A1=values["A1"] <= CRITERIA["A1"], A2=values["A2"] <= CRITERIA["A2"], A3=values["A3"] >= CRITERIA["A3"],
        A4=values["A4"] <= CRITERIA["A4"], R1=values["R1"] <= CRITERIA["R1"], R2=values["R2"] <= CRITERIA["R2"],
    )
    worst = dict(
        A2=max(((r["welch_flag_rate"], r["scenario"], r["metric"]) for r in core if r["expected"] == 0)),
        A3=min(((r["ci_coverage"], r["scenario"], r["metric"]) for r in core if np.isfinite(r["ci_coverage"]))),
        A4=max(((f["familywise_false_call_rate"], f["scenario"]) for f in fw)),
        R1=max(((r["fdr"], r["scenario"], r["metric"]) for r in with_truth)) if with_truth else None,
        R2=max(((r["fdr"], r["scenario"], r["metric"]) for r in complete_null)) if complete_null else None,
    )
    return dict(values=values, passed=passed, worst=worst)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="profileseg_validation")
    parser.add_argument("--reps", type=int, default=40)
    parser.add_argument("--templates-2cell", type=Path)
    parser.add_argument("--only", nargs="*", help="run only these designs")
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    frozen = source_hash()
    results = {}
    for name, d in designs(args.templates_2cell).items():
        if args.only and name not in args.only:
            continue
        began = time.time()
        s = d["settings"]
        config = Config(reference="WT", alternative="KO", flank=s.flank, binsize=s.binsize)
        agg = run_benchmark(reps=args.reps, settings=s, config=config, templates=d.get("templates"),
                            bootstrap=200, seed=90_000)
        reg = run_regional_benchmark(reps=args.reps, settings=s, config=config, templates=d.get("templates"), seed=95_000)
        results[name] = dict(evaluate(agg, reg), seconds=round(time.time() - began, 1),
                             aggregate=agg["summary"], regional=reg["summary"])
        print(f"{name}: {results[name]['values']}  passed={results[name]['passed']}", flush=True)
    (out / "validation.json").write_text(json.dumps(dict(source_sha256=frozen, reps=args.reps, criteria=CRITERIA,
                                                         results=results), indent=1, default=float))
    lines = ["# ProfileSeg independent validation", "", f"Frozen source SHA-256: `{frozen}`; {args.reps} runs per scenario.", "",
             "| design | A1 median false | A2 max false | A3 min coverage | A4 family-wise | R1 max FDR | R2 null FWER | R3 site power |",
             "|---|---|---|---|---|---|---|---|"]
    for name, r in results.items():
        v, p = r["values"], r["passed"]
        cell = lambda k, fmt=".3f": f"{v[k]:{fmt}}" + ("" if k == "R3" or p[k] else " ✗")  # noqa: E731
        lines.append(f"| {name} | " + " | ".join(cell(k) for k in ("A1", "A2", "A3", "A4", "R1", "R2", "R3")) + " |")
    lines += ["", "Criteria: " + ", ".join(f"{k} {'≥' if k == 'A3' else '≤'} {v}" for k, v in CRITERIA.items()), "",
              "## Worst cases", ""]
    for name, r in results.items():
        lines.append(f"* **{name}**: " + "; ".join(f"{k} {w}" for k, w in r["worst"].items()))
    (out / "VALIDATION.md").write_text("\n".join(lines) + "\n")
    print(f"report: {out / 'VALIDATION.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
