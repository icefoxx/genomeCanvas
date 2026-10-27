"""Held-out-region check of the replicate-structure explanation for yeast null-split calls.

Regions are split at random (fixed seed) into halves A and B. The label-free pairing ratio (across-pair
/ within-pair distance of per-sample meta-profiles, profileseg_replicate_structure.py) is computed on A
only; the v3 shape tests are run on B only, for every within-genotype null split (no insertion).

Prediction written before the run: splits with ratio(A) > 1.2 have higher null-call shares on B than the
others, and Spearman(ratio(A), shape_kernel null-call share on B) > 0.5.

This removes one alternative (noise in the same regions driving both the ratio and the calls). It does not
make the evidence independent of the samples: the same replicates enter both halves, so a shared
sample-level difference of any origin (batch, culture, library) is what the check can detect.

    python benchmarks/profileseg_replicate_holdout.py YEAST_CACHE.npz OUT_DIR
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from profileseg_insertion import YEAST, splits_for  # noqa: E402
from profileseg_replicate_structure import PAIRINGS, _spearman, pairing_ratios  # noqa: E402

from genomecanvas.profileseg.analysis import analyse  # noqa: E402
from genomecanvas.profileseg.data import Config, Profiles  # noqa: E402
from genomecanvas.profileseg.regional import per_region_tests  # noqa: E402

SEED = 9_600_000


def run(job):
    cache, half_b, a_ids, b_ids, name = job
    d = np.load(cache, allow_pickle=True)
    ids = list(d["sample_ids"])
    use = [ids.index(i) for i in a_ids + b_ids]
    values = d["values"][use][:, half_b]
    regions = [dict(chrom="c", start=i, end=i + 1, name=f"r{i}", strand="+") for i in range(values.shape[1])]
    smp = [dict(sample_id=i, condition="A" if k < len(a_ids) else "B") for k, i in enumerate(a_ids + b_ids)]
    reg = per_region_tests(analyse(Profiles(values, d["edges"], smp, regions, True),
                                   Config(reference="A", alternative="B", **YEAST["config"]), bootstrap=0))
    out = {m: float(np.isin(reg.tests[m]["call"], ("Higher", "Lower", "Changed")).mean()) for m in ("shape_kernel", "shape_any", "shift")}
    return name, out


def main(argv=None) -> int:
    argv = argv or sys.argv[1:]
    cache, out = argv[0], Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    git = lambda *a: subprocess.run(["git", *a], cwd=root, capture_output=True, text=True).stdout.strip()  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                regional_py_sha256=hashlib.sha256((root / "src/genomecanvas/profileseg/regional.py").read_bytes()).hexdigest()[:16],
                src_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--", "src").encode()).hexdigest()[:16],
                seed=SEED, started=datetime.datetime.now().isoformat(timespec="seconds"))
    d = np.load(cache, allow_pickle=True)
    v, ids, groups = np.nan_to_num(d["values"]), list(d["sample_ids"]), json.loads(str(d["groups"]))
    perm = np.random.default_rng(SEED).permutation(v.shape[1])
    half_a, half_b = np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2:])
    ratio = {}
    for gt, smp in groups.items():
        r = pairing_ratios(v[[ids.index(s) for s in smp]][:, half_a])
        for name in PAIRINGS:
            ratio[f"{gt}:{smp[0]}+{smp[int(name[-1]) - 1]}"] = r[name]
    jobs = [(cache, half_b, a, b, name) for name, a, b in splits_for(groups)]
    res = {}
    with ProcessPoolExecutor(6) as ex:
        for name, o in ex.map(run, jobs):
            res[name] = o
    prov["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    names = [n for n in res if n in ratio]
    r = np.array([ratio[n] for n in names])
    k = np.array([res[n]["shape_kernel"] for n in names])
    hi = r > 1.2
    summary = dict(n_splits=len(names), high_ratio_splits=int(hi.sum()),
                   null_share_kernel_high=float(k[hi].mean()) if hi.any() else float("nan"),
                   null_share_kernel_other=float(k[~hi].mean()) if (~hi).any() else float("nan"),
                   splits_with_calls_high=int((k[hi] > 0).sum()), splits_with_calls_other=int((k[~hi] > 0).sum()),
                   spearman_ratio_vs_kernel=_spearman(r, k),
                   prediction_1_higher_rate=bool(hi.any() and (~hi).any() and k[hi].mean() > k[~hi].mean()),
                   prediction_2_spearman_gt_0_5=bool(_spearman(r, k) > 0.5))
    (out / "replicate_holdout.json").write_text(json.dumps(dict(provenance=prov, regions_a=len(half_a), regions_b=len(half_b),
                                                               ratio_a=ratio, null_calls_b=res, summary=summary), indent=1))
    lines = ["# Held-out-region check: pairing ratio on half A vs null calls on half B (yeast)", "",
             f"Commit `{prov['commit'][:10]}`; regional.py `{prov['regional_py_sha256']}`; src diff `{prov['src_diff_sha256']}`; "
             f"seed {SEED}; {len(half_a)} / {len(half_b)} regions; {prov['started']} → {prov['finished']}.", "",
             "| split | pairing ratio (A) | null calls on B: shape_kernel / shape_any / shift |", "|---|---|---|"]
    for n in sorted(names, key=lambda x: -ratio[x]):
        o = res[n]
        lines.append(f"| {n} | {ratio[n]:.2f} | {100 * o['shape_kernel']:.2f} % / {100 * o['shape_any']:.2f} % / {100 * o['shift']:.2f} % |")
    lines += ["", "Pre-specified predictions: (1) higher null-call share on B for ratio(A) > 1.2: "
              f"**{'met' if summary['prediction_1_higher_rate'] else 'not met'}** "
              f"({100 * summary['null_share_kernel_high']:.2f} % vs {100 * summary['null_share_kernel_other']:.2f} %; splits with calls "
              f"{summary['splits_with_calls_high']} of {summary['high_ratio_splits']} vs {summary['splits_with_calls_other']} of "
              f"{len(names) - summary['high_ratio_splits']}); (2) Spearman > 0.5: **{'met' if summary['prediction_2_spearman_gt_0_5'] else 'not met'}** "
              f"({summary['spearman_ratio_vs_kernel']:.2f}).", ""]
    (out / "REPLICATE_HOLDOUT.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
