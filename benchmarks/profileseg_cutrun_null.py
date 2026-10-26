"""Single-cell-line 3 vs 3 real nulls and insertion on TF CUT&RUN (GSE291134, EWS::FLI1, untreated).

Two cell lines (A673; EF-dTAG clone A2.2), six replicates each of FLI1 CUT&RUN at 0 nM dTAG.
Anchors are called condition-blind from the pooled signal of all twelve samples (each scaled to
its genome-wide mean): the 5,000 strongest local maxima at 200 bp resolution, ≥ 3 kb apart, refined
to 10 bp within ±200 bp. Profiles ±1 kb in 10 bp bins; anchors are oriented along the genome (+),
so `shift` tests genomic-direction movement.

Per line: all 10 splits into 3 vs 3 (nulls; per-region and aggregate calls), and insertion of known
changes into 20 % of regions of the second group.

    python benchmarks/profileseg_cutrun_null.py --bigwigs DIR --out DIR
"""

from __future__ import annotations

import argparse
import itertools
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
from genomecanvas.profileseg.data import Config, Profiles, extract

LINES = {"A673": [f"GSM88295{49 + k}_A673_0nM_FLI1_rep{k + 1}.bigWig" for k in range(6)],
         "A2.2": [f"GSM88295{67 + k}_A2.2_0nM_FLI1_rep{k + 1}.bigWig" for k in range(6)]}
CHROMS = [f"chr{i}" for i in list(range(1, 23)) + ["X"]]
CFG = dict(anchor="center", orientation="oriented", flank=1000, binsize=10, background_width=300, core_flank=600)
INSERTIONS = ("none", "amplitude", "shift20", "shift50", "shift100", "broaden1.3")


def binned_mean(h, chrom: str, edges_bp: np.ndarray) -> np.ndarray:
    """Exact mean coverage per bin from the bigWig's intervals (cumulative area; no per-base array)."""
    iv = h.intervals(chrom)
    if not iv:
        return np.zeros(len(edges_bp) - 1)
    a = np.asarray(iv, dtype=float)
    starts, ends, vals = a[:, 0], a[:, 1], a[:, 2]
    cum = np.concatenate([[0.0], np.cumsum(vals * (ends - starts))])
    idx = np.searchsorted(ends, edges_bp, side="right")  # intervals ending at or before each edge
    area = cum[idx]
    inside = (idx < len(starts)) & (starts[np.minimum(idx, len(starts) - 1)] < edges_bp)
    j = np.minimum(idx, len(starts) - 1)
    area = area + np.where(inside, vals[j] * (edges_bp - starts[j]), 0.0)
    return np.diff(area) / np.diff(edges_bp)


def call_anchors(paths: list[Path], n: int = 5000, step: int = 200, sep: int = 3000) -> list[dict]:
    import pyBigWig
    handles = [pyBigWig.open(str(p)) for p in paths]
    scale = []
    for h in handles:  # genome-wide mean over the main chromosomes
        tot = sum((h.stats(c, type="mean", exact=False)[0] or 0) * h.chroms()[c] for c in CHROMS if c in h.chroms())
        scale.append(tot / sum(h.chroms()[c] for c in CHROMS if c in h.chroms()))
    cands = []
    for c in CHROMS:
        length = handles[0].chroms().get(c)
        if not length:
            continue
        nb = length // step
        pooled = np.zeros(nb)
        edges_bp = np.arange(nb + 1) * step
        for h, s in zip(handles, scale, strict=True):
            pooled += binned_mean(h, c, edges_bp) / s
        sm = np.convolve(pooled, np.ones(3) / 3, mode="same")
        peak = (sm[1:-1] > sm[:-2]) & (sm[1:-1] >= sm[2:])
        idx = np.flatnonzero(peak) + 1
        top = idx[np.argsort(sm[idx])[-4 * n // len(CHROMS) - 200:]]
        cands += [(float(sm[i]), c, int(i * step + step // 2)) for i in top]
    cands.sort(reverse=True)
    chosen, taken = [], {}
    for score, c, pos in cands:
        if any(abs(pos - q) < sep for q in taken.get(c, [])):
            continue
        taken.setdefault(c, []).append(pos)
        chosen.append((c, pos))
        if len(chosen) >= n:
            break
    out = []
    for c, pos in chosen:  # refine the summit at 10 bp within ±200 bp
        lo = max(pos - 200, 0)
        prof = sum(np.nan_to_num(np.array(h.values(c, lo, pos + 200), dtype=float)).reshape(-1, 10).mean(1) / s
                   for h, s in zip(handles, scale, strict=True))
        summit = lo + int(np.argmax(prof)) * 10 + 5
        out.append(dict(chrom=c, start=summit, end=summit + 1, name=f"{c}:{summit}", strand="+"))
    return out


def run(job):
    cache, line, a_ids, b_ids, insertion, rep = job
    d = np.load(cache, allow_pickle=True)
    ids = list(d["sample_ids"])
    use = [ids.index(i) for i in a_ids + b_ids]
    values, edges = d["values"][use], d["edges"]
    cfg = Config(reference="A", alternative="B", **CFG)
    rng = np.random.default_rng(zlib.crc32(f"{line}|{a_ids}|{insertion}|{rep}".encode()))
    bg = M.background(values, edges, cfg.background_width, "mean")
    core = M.core_mask(edges, cfg.core_flank)
    cand = np.flatnonzero(np.nan_to_num((values - bg[..., None])[..., core]).sum(-1).mean(0) > 0)
    chosen = rng.choice(cand, int(0.2 * len(cand)), replace=False) if insertion != "none" else np.array([], int)
    if insertion != "none":
        values = I.insert(values, edges, list(range(3, 6)), chosen, insertion, rng, cfg.background_width)
    truth = np.zeros(values.shape[1], bool)
    truth[chosen] = insertion not in ("none", "amplitude")
    smp = [dict(sample_id=i, condition="A" if i in a_ids else "B") for i in a_ids + b_ids]
    regions = [dict(chrom="c", start=i, end=i + 1, name=f"r{i}", strand="+") for i in range(values.shape[1])]
    a = analyse(Profiles(values, edges, smp, regions, True), cfg, bootstrap=0)
    from genomecanvas.profileseg.regional import per_region_tests
    reg = per_region_tests(a)
    out = {}
    for m in ("shape_kernel", "shift"):
        t = reg.tests[m]
        called = np.isin(t["call"], ("Higher", "Lower", "Changed"))
        tr = truth[t["region_index"]]
        nc = int(called.sum())
        out[m] = dict(called=float(called.mean()), any=float(called.any()),
                      fdp=float((called & ~tr).sum() / nc) if nc else 0.0,
                      power=float((called & tr).sum() / tr.sum()) if tr.any() else float("nan"))
    # aggregate report: primary and Holm-adjusted secondary calls (comparators excluded)
    calls = [e for e in a.effects if e.get("role") in ("primary", "secondary") and e["metric"] != "group_w1"]
    out["aggregate_any_call"] = float(any(e["call"] in ("higher", "lower") for e in calls))
    out["aggregate_primary_call"] = float(any(e["call"] in ("higher", "lower") for e in calls if e.get("role") == "primary"))
    return line, "+".join(a_ids), insertion, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bigwigs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=1)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    bw = Path(args.bigwigs)
    cache = out / "cutrun_profiles.npz"
    if not cache.exists():
        paths = [bw / f for fs in LINES.values() for f in fs]
        regions = call_anchors(paths)
        samples = [dict(sample_id=f"{line}_{k + 1}", condition=line, signal_file=str(bw / f))
                   for line, fs in LINES.items() for k, f in enumerate(fs)]
        prof = extract(Config(**CFG), samples, regions)
        np.savez_compressed(cache, values=prof.values, edges=prof.edges, sample_ids=np.array([s["sample_id"] for s in samples]),
                            anchors=np.array([r["name"] for r in regions]))
    jobs = []
    for line in LINES:
        ids = [f"{line}_{k + 1}" for k in range(6)]
        for a in itertools.combinations(ids, 3):
            if ids[0] not in a:
                continue
            b = [x for x in ids if x not in a]
            jobs += [(str(cache), line, list(a), b, ins, r) for ins in INSERTIONS for r in range(args.reps)]
    rows = []
    with ProcessPoolExecutor(6) as ex:
        for line, split, ins, res in ex.map(run, jobs):
            rows.append(dict(line=line, split=split, insertion=ins, **res))
    (out / "cutrun.json").write_text(json.dumps(rows, indent=1))
    lines = ["# TF CUT&RUN (GSE291134): single-line 3 vs 3 nulls and insertion", ""]
    for line in LINES:
        lines += [f"## {line} (10 splits)", "",
                  "| test | null: share called / FWER | amplitude-only | " + " | ".join(f"{i}: power / FDP" for i in INSERTIONS[2:]) + " |",
                  "|---|---|---|" + "---|" * (len(INSERTIONS) - 2)]
        for m in ("shape_kernel", "shift"):
            sel = lambda ins: [r[m] for r in rows if r["line"] == line and r["insertion"] == ins]
            mean = lambda ins, k: float(np.nanmean([x[k] for x in sel(ins)]))
            lines.append(f"| {m} | {mean('none', 'called'):.4f} / {mean('none', 'any'):.2f} | {mean('amplitude', 'called'):.4f} | "
                         + " | ".join(f"{mean(i, 'power'):.2f} / {mean(i, 'fdp'):.3f}" for i in INSERTIONS[2:]) + " |")
        nulls = [r for r in rows if r["line"] == line and r["insertion"] == "none"]
        lines += ["", f"Aggregate report under the null: any primary/secondary call in "
                  f"{np.mean([r['aggregate_any_call'] for r in nulls]):.2f} of splits; primary calls in "
                  f"{np.mean([r['aggregate_primary_call'] for r in nulls]):.2f}.", "",
                  "Per-split null calls (shape_kernel / shift share of regions):", ""]
        lines += [f"* {r['split']}: {r['shape_kernel']['called']:.4f} / {r['shift']['called']:.4f}" for r in nulls]
        lines.append("")
    (out / "CUTRUN.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
