"""In-silico insertion on real replicates: per-region FDR / power with real noise and real shapes.

Base experiments are splits of real replicates of one condition (no true difference). Known changes
are inserted into a random 20 % of regions, in the second group only: shift (20 / 50 / 100 bp,
random sign per region), broadening (×1.3 / ×1.5 about the anchor, mass preserved) or an
amplitude-only change (background-corrected signal ×2; a negative control for shape tests). The
noise moves with the signal, as it would for a real change.

Methods compared per region (BH 5 %):
  ProfileSeg `shape_kernel`, `shift`                    (current)
  v2 centroid, moderated t                              (previous ProfileSeg version)
  Welch t on centroid, no variance shrinkage            (common replicate-level practice)
  pausing index (log promoter / gene-body signal), Welch (PRO-seq practice; shape and amplitude mixed)
and one global test: Wilcoxon signed-rank over regions of the group-mean centroid difference (regions
as units, the pseudo-replicated practice).

Input caches are built by --build-yeast / --build-proseq from GEO bigWigs (read with pyBigWig).

    python benchmarks/profileseg_insertion.py --out DIR --yeast-cache Y.npz --proseq-cache P.npz
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
from scipy import stats

from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, Profiles
from genomecanvas.profileseg.regional import bh, moderated_t, per_region_tests

INSERTIONS = ("none", "amplitude", "shift20", "shift50", "shift100", "broaden1.3", "broaden1.5")
FRACTION = 0.2

YEAST = dict(config=dict(anchor="start", orientation="oriented", flank=1500, binsize=10, background_width=300,
                         core_flank=1000))
PROSEQ = dict(config=dict(anchor="start", orientation="oriented", flank=1000, binsize=10, background_width=300,
                          core_flank=500))


def splits_for(groups: dict, balanced_by: dict | None = None) -> list[tuple[str, list, list]]:
    """2 vs 2 splits within each 4-replicate group; for PRO-seq, 3 vs 3 splits with one replicate of
    each cell line per group (line differences cancel)."""
    out = []
    if balanced_by is None:
        for g, ids in groups.items():
            for a in itertools.combinations(ids, len(ids) // 2):
                if ids[0] not in a:
                    continue
                out.append((f"{g}:{'+'.join(a)}", list(a), [i for i in ids if i not in a]))
    else:
        lines = list(balanced_by)
        for pick in itertools.product(*[range(2) for _ in lines]):
            if pick[0] != 0:
                continue
            a = [balanced_by[l][k] for l, k in zip(lines, pick, strict=True)]
            b = [balanced_by[l][1 - k] for l, k in zip(lines, pick, strict=True)]
            out.append(("+".join(a), a, b))
    return out


def insert(values: np.ndarray, edges: np.ndarray, alt: list[int], regions: np.ndarray, kind: str, rng,
           background_width: int) -> np.ndarray:
    v = values.copy()
    c = 0.5 * (edges[:-1] + edges[1:])
    nb = max(1, int(background_width / float(np.median(np.diff(edges)))))
    for r in regions:
        sign = rng.choice([-1, 1])
        for s in alt:
            x = np.nan_to_num(v[s, r])
            bg = np.r_[x[:nb], x[-nb:]].mean()
            if kind.startswith("shift"):
                d = float(kind[5:]) * sign
                v[s, r] = np.interp(c - d, c, x, left=x[0], right=x[-1])
            elif kind.startswith("broaden"):
                f = float(kind[7:])
                y = bg + np.interp(c / f, c, x - bg)
                v[s, r] = bg + (y - bg) * max((x - bg).sum(), 1e-9) / max((y - bg).sum(), 1e-9)
            elif kind == "amplitude":
                v[s, r] = bg + (x - bg) * 2.0
    return v


def pausing_index(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    c = 0.5 * (edges[:-1] + edges[1:])
    prom = np.nansum(values[..., (c >= 0) & (c < 250)], -1)
    body = np.nansum(values[..., (c >= 400) & (c < 1000)], -1)
    return np.log2((prom + 1) / (body + 1))


def run(job):
    cache, kind_of_data, name, a_ids, b_ids, insertion, rep = job
    d = np.load(cache, allow_pickle=True)
    ids = list(d["sample_ids"])
    use = [ids.index(i) for i in a_ids + b_ids]
    values, edges = d["values"][use], d["edges"]
    cfgd = (YEAST if kind_of_data == "yeast" else PROSEQ)["config"]
    cfg = Config(reference="A", alternative="B", **cfgd)
    regions = [dict(chrom="c", start=i, end=i + 1, name=f"r{i}", strand="+") for i in range(values.shape[1])]
    rng = np.random.default_rng(zlib.crc32(f"{name}|{insertion}|{rep}".encode()))  # deterministic across runs
    n_a = len(a_ids)
    alt = list(range(n_a, len(use)))
    # candidate regions: positive pooled background-corrected core signal (the shape-testable set)
    bg = M.background(values, edges, cfg.background_width, "mean")
    core = M.core_mask(edges, cfg.core_flank)
    pooled = np.nan_to_num((values - bg[..., None])[..., core]).sum(-1).mean(0)
    cand = np.flatnonzero(pooled > 0)
    chosen = rng.choice(cand, int(FRACTION * len(cand)), replace=False) if insertion != "none" else np.array([], int)
    if insertion != "none":
        values = insert(values, edges, alt, chosen, insertion, rng, cfg.background_width)
    truth = np.zeros(values.shape[1], bool)
    if insertion not in ("none", "amplitude"):
        truth[chosen] = True
    smp = [dict(sample_id=i, condition="A" if k < n_a else "B") for k, i in enumerate(a_ids + b_ids)]
    a = analyse(Profiles(values, edges, smp, regions, True), cfg, bootstrap=0)
    reg = per_region_tests(a)
    out = {}

    def score(idx, called):
        t = truth[idx]
        nc = int(called.sum())
        return dict(called=float(called.mean()), any=float(called.any()),
                    fdp=float((called & ~t).sum() / nc) if nc else 0.0,
                    power=float((called & t).sum() / t.sum()) if t.any() else float("nan"))

    for m in ("shape_kernel", "shift", "centroid"):
        if m in reg.tests:
            t = reg.tests[m]
            out["v2 centroid (moderated t)" if m == "centroid" else m] = score(
                t["region_index"], np.isin(t["call"], ("Higher", "Lower", "Changed")))
    # Welch on centroid per region, no shrinkage
    idx = reg.tests["centroid"]["region_index"]
    p_full = M.normalise(np.clip(np.nan_to_num(a.core_signal[:, idx]), 0, None) * M.widths(a.core_edges))
    cen = M.centroid(p_full, a.core_edges)
    ref, altx = np.arange(n_a), np.arange(n_a, len(use))
    w = stats.ttest_ind(cen[altx], cen[ref], axis=0, equal_var=False)
    out["Welch t on centroid"] = score(idx, bh(np.nan_to_num(w.pvalue, nan=1.0)) < 0.05)
    # global, regions as units
    diff = cen[altx].mean(0) - cen[ref].mean(0)
    out["global Wilcoxon over regions"] = dict(called=float(stats.wilcoxon(diff).pvalue < 0.05))
    if kind_of_data == "proseq":
        pi = pausing_index(values, edges)
        ok = np.isfinite(pi).all(0)
        w = stats.ttest_ind(pi[altx][:, ok], pi[ref][:, ok], axis=0, equal_var=False)
        out["pausing index (Welch)"] = score(np.flatnonzero(ok), bh(np.nan_to_num(w.pvalue, nan=1.0)) < 0.05)
        mt = moderated_t(pi[:, ok], ref, altx)
        out["pausing index (moderated t)"] = score(np.flatnonzero(ok), bh(mt["p"]) < 0.05)
    return kind_of_data, name, insertion, out


def build_yeast(directory: Path, bed: Path, out: Path) -> None:
    from genomecanvas.profileseg.data import extract, read_regions
    groups = {"WT": ("wildtype", range(8603534, 8603538)), "dKO": ("isw1_null_chd1_null", range(8603568, 8603572)),
              "rtf1-8-9A": ("rtf1_8_9_A", range(8603538, 8603542)), "rtf1-11A": ("rtf1_11_A", range(8603542, 8603546)),
              "rtf1": ("rtf1_null", range(8603549, 8603553)), "isw1rtf1": ("isw1_null_rtf1_null", range(8603564, 8603568))}
    samples, gmap = [], {}
    for g, (tag, gsms) in groups.items():
        for k, gsm in enumerate(gsms, 1):
            sid = f"{g}_{k}"
            samples.append(dict(sample_id=sid, condition=g,
                                signal_file=str(directory / f"GSM{gsm}_{tag}_bio{k}_mononucleosomes_CPM.bw")))
            gmap.setdefault(g, []).append(sid)
    cfg = Config(**YEAST["config"])
    prof = extract(cfg, samples, read_regions(bed))
    np.savez_compressed(out, values=prof.values, edges=prof.edges, sample_ids=np.array([s["sample_id"] for s in samples]),
                        groups=json.dumps(gmap))


def build_proseq(directory: Path, out: Path) -> None:
    import gzip
    from genomecanvas.profileseg.data import extract
    files = {"NELFB_1": "GSM7472402_NELFB-DMSO-{}-dedup-Rep1.bw", "NELFB_2": "GSM7472403_NELFB-DMSO-{}-dedup-Rep2.bw",
             "PAF1_1": "GSM7472394_PAF1-DMSO-{}-dedup-Rep1.bw", "PAF1_2": "GSM7472395_PAF1-DMSO-{}-dedup-Rep2.bw",
             "RTF1_1": "GSM7472405_RTF1-DMSO-{}-dedup-Rep1.bw", "RTF1_2": "GSM7472406_RTF1-DMSO-{}-dedup-Rep2.bw"}
    rows = [l.split("\t") for l in gzip.open(directory / "hg38_ncbiRefSeqSelect.txt.gz", "rt")]
    main = {f"chr{i}" for i in list(range(1, 23)) + ["X"]}
    regions, seen = {"+": [], "-": []}, set()
    for r in rows:
        if not r[1].startswith("NM_") or r[2] not in main or int(r[5]) - int(r[4]) < 2000 or r[12] in seen:
            continue
        seen.add(r[12])
        tss = int(r[4]) if r[3] == "+" else int(r[5]) - 1
        regions[r[3]].append(dict(chrom=r[2], start=tss, end=tss + 1, name=r[12], strand=r[3]))
    cfg = Config(**PROSEQ["config"])
    blocks = []
    for strand, kind in (("+", "FW"), ("-", "RV")):
        smp = [dict(sample_id=k, condition=k.split("_")[0], signal_file=str(directory / f.format(kind))) for k, f in files.items()]
        p = extract(cfg, smp, regions[strand])
        blocks.append(p.values)
    values = np.concatenate(blocks, axis=1)
    # genes with promoter-proximal signal in the pooled data (expressed), condition-blind
    c = 0.5 * (p.edges[:-1] + p.edges[1:])
    prom = np.nansum(values[..., (c >= 0) & (c < 300)], -1).mean(0)
    keep = prom > np.quantile(prom, 0.5)
    np.savez_compressed(out, values=values[:, keep], edges=p.edges, sample_ids=np.array(list(files)),
                        groups=json.dumps({"NELFB": ["NELFB_1", "NELFB_2"], "PAF1": ["PAF1_1", "PAF1_2"],
                                           "RTF1": ["RTF1_1", "RTF1_2"]}))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--yeast-cache")
    ap.add_argument("--proseq-cache")
    ap.add_argument("--build-yeast", nargs=2, metavar=("BIGWIG_DIR", "BED"))
    ap.add_argument("--build-proseq", metavar="BIGWIG_DIR")
    ap.add_argument("--reps", type=int, default=2, help="random insertion sets per base split")
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    if args.build_yeast:
        build_yeast(Path(args.build_yeast[0]), Path(args.build_yeast[1]), Path(args.yeast_cache))
    if args.build_proseq:
        build_proseq(Path(args.build_proseq), Path(args.proseq_cache))
    jobs = []
    if args.yeast_cache:
        g = json.loads(str(np.load(args.yeast_cache, allow_pickle=True)["groups"]))
        for name, a, b in splits_for(g):
            jobs += [(args.yeast_cache, "yeast", name, a, b, ins, r) for ins in INSERTIONS for r in range(args.reps)]
    if args.proseq_cache:
        g = json.loads(str(np.load(args.proseq_cache, allow_pickle=True)["groups"]))
        for name, a, b in splits_for(g, balanced_by=g):
            jobs += [(args.proseq_cache, "proseq", name, a, b, ins, r) for ins in INSERTIONS for r in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for data, name, ins, res in pool.map(run, jobs):
            for m, v in res.items():
                agg.setdefault(data, {}).setdefault(m, {}).setdefault(ins, []).append(v)
    lines = ["# In-silico insertion on real replicates", ""]
    table = {}
    for data, methods in agg.items():
        lines += [f"## {data}", "", "| method | null: share called / FWER | amplitude-only: share called | "
                  + " | ".join(f"{i}: power / FDP" for i in INSERTIONS[2:]) + " |", "|---|---|---|" + "---|" * (len(INSERTIONS) - 2)]
        for m, by in methods.items():
            mean = lambda ins, k: float(np.nanmean([r[k] for r in by[ins]])) if k in by[ins][0] else float("nan")
            table.setdefault(data, {})[m] = {ins: {k: mean(ins, k) for k in by[ins][0]} for ins in by}
            if m.startswith("global"):
                cells = " | ".join(f"{mean(i, 'called'):.2f} (share of runs significant)" for i in INSERTIONS[2:])
                lines.append(f"| {m} | {mean('none', 'called'):.2f} of runs significant | {mean('amplitude', 'called'):.2f} | {cells} |")
                continue
            cells = " | ".join(f"{mean(i, 'power'):.2f} / {mean(i, 'fdp'):.3f}" for i in INSERTIONS[2:])
            lines.append(f"| {m} | {mean('none', 'called'):.4f} / {mean('none', 'any'):.2f} | {mean('amplitude', 'called'):.4f} | {cells} |")
        lines.append("")
    (out / "insertion.json").write_text(json.dumps(table, indent=1))
    (out / "INSERTION.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
