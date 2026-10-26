"""ProfileSeg vs MMDiff2 on the same real BAMs (yeast MNase-seq, GSE280671; reads from ENA).

BAMs: the first 3 M read pairs per sample, bowtie2 to sacCer3, MAPQ ≥ 10, proper pairs. Both tools
see the same fragments: ProfileSeg gets fragment-midpoint (dyad) counts in 10 bp bins over
ATG ± 1.5 kb, extracted with samtools; MMDiff2 reads the first mates in single-end mode (its
paired-end mode fails in version 1.30.0, see MMD_R).

Comparisons (2 vs 2): real nulls WT 1+2 vs 3+4 and WT 1+3 vs 2+4; positive control WT 1+2 vs
isw1Δ chd1Δ 1+2. MMDiff2 runs on a fixed random subset of genes (it is slow); ProfileSeg on the same
subset and on all genes.

    python benchmarks/profileseg_realbam.py --bams DIR --bed orf_atg.bed --out DIR [--genes 300]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, Profiles, read_regions
from genomecanvas.profileseg.regional import per_region_tests

# samtools: $SAMTOOLS, else on PATH, else the conda env used for these benchmarks
SAMTOOLS = Path(os.environ.get("SAMTOOLS") or shutil.which("samtools") or Path.home() / "anaconda3/envs/psbench/bin/samtools")
FLANK, BIN, CORE = 1500, 10, 1000
COMPARISONS = {"null WT 1+2 vs 3+4": (["WT_rep1", "WT_rep2"], ["WT_rep3", "WT_rep4"]),
               "null WT 1+3 vs 2+4": (["WT_rep1", "WT_rep3"], ["WT_rep2", "WT_rep4"]),
               "WT 1+2 vs dKO 1+2": (["WT_rep1", "WT_rep2"], ["dKO_rep1", "dKO_rep2"])}

MMD_R = r"""
suppressMessages({library(MMDiff2); library(GenomicRanges)})
args <- commandArgs(TRUE); run <- args[1]
peaks <- read.table(file.path(run, "peaks.tsv"), header = TRUE)
gr <- GRanges(peaks$chrom, IRanges(peaks$start + 1, peaks$end)); names(gr) <- peaks$name
meta <- list(ExpData = list(dataDir = run, sampleSheet = "samples.csv"))
md <- DBAmmd(meta); md <- setRegions(md, gr)
# MMDiff2 1.30.0's paired-end branch of getReads references an undefined RawCounts.p and stops,
# so first mates are analysed in single-end mode (fragment centres estimated by MMDiff2).
md <- suppressMessages(getPeakReads(md, pairedEnd = FALSE, run.parallel = FALSE))
md <- suppressMessages(estimateFragmentCenters(md))
md <- compHists(md, bin.length = 20, whichPos = "Center")
md <- suppressMessages(compDists(md, dist.method = "MMD", run.parallel = FALSE))
md <- setContrast(md, contrast = "byCondition")
md <- suppressMessages(compPvals(md, dist.method = "MMD"))
tab <- md@Contrasts[[1]]$MMD.locfit$mmd.Table
writeLines(jsonlite::toJSON(list(p = as.numeric(tab[, "pval"]), padj = as.numeric(tab[, "padj"])), digits = NA),
           file.path(run, "mmd.json"))
"""


def midpoint_profiles(bams: dict, regions: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    """Fragment-midpoint counts per 10 bp bin over anchor ± FLANK (strand-oriented)."""
    edges = np.arange(-FLANK, FLANK + BIN, BIN, dtype=float)
    out = np.zeros((len(bams), len(regions), len(edges) - 1))
    for i, (sid, bam) in enumerate(bams.items()):
        bed = Path(bam).with_suffix(".regions.bed")
        bed.write_text("".join(f"{r['chrom']}\t{max(r['start'] - FLANK - 600, 0)}\t{r['start'] + FLANK + 600}\n" for r in regions))
        # one pass over the needed regions: first mate of each proper pair, with TLEN
        txt = subprocess.run([str(SAMTOOLS), "view", "-f", "66", "-L", str(bed), bam], capture_output=True, text=True).stdout
        by_chrom: dict = {}
        for line in txt.splitlines():
            f = line.split("\t", 9)
            pos, pnext, tlen = int(f[3]), int(f[7]), abs(int(f[8]))
            mid = min(pos, pnext) + tlen // 2  # leftmost mate start + half the fragment
            by_chrom.setdefault(f[2], []).append(mid)
        mids = {c: np.sort(np.array(v)) for c, v in by_chrom.items()}
        for k, r in enumerate(regions):
            m = mids.get(r["chrom"])
            if m is None:
                continue
            lo, hi = np.searchsorted(m, [r["start"] - FLANK, r["start"] + FLANK])
            rel = m[lo:hi] - r["start"]
            if r["strand"] == "-":
                rel = -rel
            out[i, k] = np.histogram(rel, bins=edges)[0]
    return out, edges


def run_profileseg(values, edges, regions, ids, a_ids, b_ids):
    use = [ids.index(i) for i in a_ids + b_ids]
    smp = [dict(sample_id=i, condition="A" if i in a_ids else "B") for i in a_ids + b_ids]
    cfg = Config(reference="A", alternative="B", anchor="start", orientation="oriented", flank=FLANK, binsize=BIN,
                 background_width=300, core_flank=CORE)
    t = time.perf_counter()
    a = analyse(Profiles(values[use], edges, smp, regions, True), cfg, bootstrap=0)
    reg = per_region_tests(a)
    dt = time.perf_counter() - t
    out = {"seconds": dt}
    for m in ("shape_kernel", "shift"):
        c = np.isin(reg.tests[m]["call"], ("Higher", "Lower", "Changed"))
        out[m] = (int(c.sum()), int(len(c)))
    if "shift" in reg.tests:
        calls = reg.tests["shift"]["call"]
        out["shift_direction"] = (int((calls == "Higher").sum()), int((calls == "Lower").sum()))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--genes", type=int, default=300)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    bams = {f"{g}_rep{k}": str(Path(args.bams) / f"{g}_rep{k}.bam") for g in ("WT", "dKO") for k in range(1, 5)}
    regions = read_regions(Path(args.bed))
    cache = out / "midpoints.npz"
    if not cache.exists():
        values, edges = midpoint_profiles(bams, regions)
        np.savez_compressed(cache, values=values, edges=edges)
    d = np.load(cache); values, edges = d["values"], d["edges"]
    ids = list(bams)
    # fixed random subset among genes with signal in every sample (condition-blind)
    core = np.abs(0.5 * (edges[:-1] + edges[1:])) < CORE
    depth = values[..., core].sum(-1)
    eligible = np.flatnonzero((depth > 50).all(0))
    subset = np.sort(np.random.default_rng(11).choice(eligible, min(args.genes, len(eligible)), replace=False))
    results = {}
    script = out / "mmd.R"; script.write_text(MMD_R)
    for label, (a_ids, b_ids) in COMPARISONS.items():
        res = {"ProfileSeg all genes": run_profileseg(values, edges, regions, ids, a_ids, b_ids),
               "ProfileSeg subset": run_profileseg(values[:, subset], edges, [regions[i] for i in subset], ids, a_ids, b_ids)}
        run_dir = out / label.replace(" ", "_").replace("+", "p")
        run_dir.mkdir(exist_ok=True)
        rows = ["SampleID,Condition,Replicate,bamReads"]
        for k, sid in enumerate(a_ids + b_ids):
            r1 = Path(bams[sid]).with_suffix(".R1.bam")
            if not r1.exists():
                subprocess.run([str(SAMTOOLS), "view", "-b", "-f", "64", "-o", str(r1), bams[sid]], check=True)
                subprocess.run([str(SAMTOOLS), "index", str(r1)], check=True)
            rows.append(f"{sid},{1 if sid in a_ids else 2},{k % 2 + 1},{r1}")
        (run_dir / "samples.csv").write_text("\n".join(rows) + "\n")
        (run_dir / "peaks.tsv").write_text("chrom\tstart\tend\tname\n" + "".join(
            f"{regions[i]['chrom']}\t{regions[i]['start'] - CORE}\t{regions[i]['start'] + CORE}\t{regions[i]['name']}\n" for i in subset))
        t = time.perf_counter()
        proc = subprocess.run(["/usr/bin/time", "-l", "Rscript", str(script), str(run_dir)], capture_output=True, text=True)
        dt = time.perf_counter() - t
        if proc.returncode:
            res["MMDiff2"] = dict(error=proc.stderr[-800:])
        else:
            mm = json.loads((run_dir / "mmd.json").read_text())
            p, padj = np.array(mm["p"], float), np.array(mm["padj"], float)
            from genomecanvas.profileseg.regional import bh
            res["MMDiff2"] = dict(seconds=dt, bh05=int(np.sum(np.nan_to_num(bh(p), nan=1) < 0.05)),
                                  own_padj01=int(np.sum(np.nan_to_num(padj, nan=1) < 0.1)), tested=int(np.isfinite(p).sum()),
                                  raw_p05=float(np.nanmean(p < 0.05)))
        results[label] = res
        print(label, json.dumps(res), flush=True)
    (out / "realbam.json").write_text(json.dumps(results, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
