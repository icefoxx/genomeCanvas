"""Batch-aware label assignments on real yeast MNase replicates (GSE280671).

GEO does not record processing batches; batches are taken as replicate pairs {1,2} and {3,4},
the structure found in benchmarks/profileseg_insertion.py (cross-pair / within-pair meta-profile
distance 1.5–2.7 in three genotypes). For each genotype: the batch-confounded split (1+2 vs 3+4)
and the two batch-balanced splits, with batch handling off, default (refuse confounded designs) and
within (permute within batches); null and insertions (20 % of regions, second group).

    python benchmarks/profileseg_batch_aware.py --cache ins_yeast.npz --out DIR
"""

from __future__ import annotations

import argparse
import json
import os
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

import profileseg_insertion as I
from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, Profiles
from genomecanvas.profileseg.regional import per_region_tests

INSERTIONS = ("none", "shift20", "shift50", "broaden1.3")


def run(job):
    cache, genotype, a_ids, b_ids, insertion, rep, aware = job
    d = np.load(cache, allow_pickle=True)
    ids = list(d["sample_ids"])
    use = [ids.index(i) for i in a_ids + b_ids]
    values, edges = d["values"][use], d["edges"]
    mode = {"off": dict(batch_aware=False), "default": dict(batch_aware=True),
            "within": dict(batch_aware=True, batch_permutation="within")}[aware]
    cfg = Config(reference="A", alternative="B", **mode, **I.YEAST["config"])
    rng = np.random.default_rng(zlib.crc32(f"{genotype}|{a_ids}|{insertion}|{rep}".encode()))
    bg = M.background(values, edges, cfg.background_width, "mean")
    core = M.core_mask(edges, cfg.core_flank)
    cand = np.flatnonzero(np.nan_to_num((values - bg[..., None])[..., core]).sum(-1).mean(0) > 0)
    chosen = rng.choice(cand, int(0.2 * len(cand)), replace=False) if insertion != "none" else np.array([], int)
    if insertion != "none":
        values = I.insert(values, edges, list(range(2, 4)), chosen, insertion, rng, cfg.background_width)
    truth = np.zeros(values.shape[1], bool); truth[chosen] = insertion != "none"
    batch = lambda sid: "b1" if sid.endswith(("_1", "_2")) else "b2"
    smp = [dict(sample_id=i, condition="A" if i in a_ids else "B", batch=batch(i)) for i in a_ids + b_ids]
    regions = [dict(chrom="c", start=i, end=i + 1, name=f"r{i}", strand="+") for i in range(values.shape[1])]
    reg = per_region_tests(analyse(Profiles(values, edges, smp, regions, True), cfg, bootstrap=0))
    out = {}
    for m in ("shape_kernel", "shift"):
        t = reg.tests[m]
        tested = t["call"] != "not_tested"
        called = np.isin(t["call"], ("Higher", "Lower", "Changed"))
        tr = truth[t["region_index"]]
        nc = int(called.sum())
        out[m] = dict(tested=float(tested.mean()), called=float(called.mean()),
                      fdp=float((called & ~tr).sum() / nc) if nc else 0.0,
                      power=float((called & tr).sum() / tr.sum()) if tr.any() else float("nan"))
    return genotype, "+".join(a_ids), insertion, aware, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=2)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    groups = json.loads(str(np.load(args.cache, allow_pickle=True)["groups"]))
    jobs = []
    for g, s in groups.items():
        for a in ([s[0], s[1]], [s[0], s[2]], [s[0], s[3]]):
            b = [x for x in s if x not in a]
            for ins in INSERTIONS:
                for aware in ("off", "default", "within"):
                    jobs += [(args.cache, g, a, b, ins, r, aware) for r in range(args.reps if ins != "none" else 1)]
    rows = []
    with ProcessPoolExecutor(5) as ex:
        for g, split, ins, aware, res in ex.map(run, jobs):
            rows.append(dict(genotype=g, split=split, confounded=split.endswith("_2"), insertion=ins, aware=aware, **res))
    (out / "batch_aware.json").write_text(json.dumps(rows, indent=1))
    lines = ["# Batch-aware label assignments (yeast MNase, 2 vs 2)", ""]
    for conf in (True, False):
        lines += [f"## {'batch-confounded split (1+2 vs 3+4)' if conf else 'batch-balanced splits'}", "",
                  "| batch mode | test | share tested | null: share called | " + " | ".join(f"{i}: power / FDP" for i in INSERTIONS[1:]) + " |",
                  "|---|---|---|---|" + "---|" * (len(INSERTIONS) - 1)]
        for aware in ("off", "default", "within"):
            for m in ("shape_kernel", "shift"):
                sel = lambda ins: [r[m] for r in rows if r["confounded"] == conf and r["aware"] == aware and r["insertion"] == ins]
                mean = lambda ins, k: float(np.nanmean([x[k] for x in sel(ins)]))
                cells = " | ".join(f"{mean(i, 'power'):.2f} / {mean(i, 'fdp'):.3f}" for i in INSERTIONS[1:])
                lines.append(f"| {aware} | {m} | {mean('none', 'tested'):.2f} | {mean('none', 'called'):.4f} | {cells} |")
        lines.append("")
    (out / "BATCH_AWARE.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
