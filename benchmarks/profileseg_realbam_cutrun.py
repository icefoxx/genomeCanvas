"""ProfileSeg vs MMDiff2 on real TF CUT&RUN BAMs (GSE291134, A673, six untreated replicates).

Reads: first 4 M pairs per replicate from ENA, bowtie2 to hg38 chr1 + chr2 only (MAPQ ≥ 30 removes
reads from other chromosomes), proper pairs. Anchors: the condition-blind FLI1 summits of
benchmarks/profileseg_cutrun_null.py that fall on chr1/chr2.

Comparisons (3 vs 3): real nulls {1,2,3} vs {4,5,6} and {1,3,5} vs {2,4,6}; a read-level insertion on
the second split: in replicates 2, 4 and 6, every fragment whose midpoint lies within ±1 kb of a
chosen anchor (20 % of the tested anchors) is moved by +60 bp. Both tools see the same fragments:
ProfileSeg as fragment-midpoint counts (10 bp bins, ±1 kb); MMDiff2 1.30.0 as first mates in
single-end mode (its paired-end mode fails in this version).

    python benchmarks/profileseg_realbam_cutrun.py --bams DIR --anchors cutrun_profiles.npz --out DIR
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, Profiles
from genomecanvas.profileseg.regional import bh, per_region_tests
from profileseg_realbam import MMD_R, SAMTOOLS

FLANK, BIN, CORE, SHIFT = 1000, 10, 600, 60
SPLITS = {"null 1,2,3 vs 4,5,6": ([1, 2, 3], [4, 5, 6]), "null 1,3,5 vs 2,4,6": ([1, 3, 5], [2, 4, 6]),
          "insertion +60 bp, 1,3,5 vs 2,4,6": ([1, 3, 5], [2, 4, 6])}


def read_pairs(bam: str, bed: Path) -> list[tuple]:
    """(chrom, leftmost start, fragment length) of each proper pair (first mates)."""
    txt = subprocess.run([str(SAMTOOLS), "view", "-f", "66", "-L", str(bed), bam], capture_output=True, text=True).stdout
    out = []
    for line in txt.splitlines():
        f = line.split("\t", 9)
        out.append((f[2], min(int(f[3]), int(f[7])), abs(int(f[8]))))
    return out


def midpoints(bams: dict, anchors: list[tuple], bed: Path) -> np.ndarray:
    edges = np.arange(-FLANK, FLANK + BIN, BIN, dtype=float)
    vals = np.zeros((len(bams), len(anchors), len(edges) - 1))
    for i, bam in enumerate(bams.values()):
        mids: dict = {}
        for c, left, tlen in read_pairs(bam, bed):
            mids.setdefault(c, []).append(left + tlen // 2)
        mids = {c: np.sort(np.array(v)) for c, v in mids.items()}
        for k, (c, pos) in enumerate(anchors):
            m = mids.get(c, np.array([]))
            lo, hi = np.searchsorted(m, [pos - FLANK, pos + FLANK])
            vals[i, k] = np.histogram(m[lo:hi] - pos, bins=edges)[0]
    return vals


def shifted_bam(bam: str, chosen: list[tuple], out: Path) -> str:
    """Copy of `bam` with fragments near chosen anchors moved by SHIFT bp (both mates)."""
    by = {}
    for c, pos in chosen:
        by.setdefault(c, []).append(pos)
    by = {c: np.sort(np.array(v)) for c, v in by.items()}
    sam = subprocess.run([str(SAMTOOLS), "view", "-h", bam], capture_output=True, text=True).stdout.splitlines()
    lines = []
    for line in sam:
        if line.startswith("@"):
            lines.append(line)
            continue
        f = line.split("\t")
        pos, pnext, tlen = int(f[3]), int(f[7]), abs(int(f[8]))
        mid = min(pos, pnext) + tlen // 2
        a = by.get(f[2])
        if a is not None and tlen > 0:
            j = np.searchsorted(a, mid)
            near = (j < len(a) and abs(a[j] - mid) <= FLANK) or (j > 0 and abs(a[j - 1] - mid) <= FLANK)
            if near:
                f[3], f[7] = str(pos + SHIFT), str(pnext + SHIFT)
        lines.append("\t".join(f))
    tmp = out.with_suffix(".sam")
    tmp.write_text("\n".join(lines) + "\n")
    subprocess.run([str(SAMTOOLS), "sort", "-o", str(out), str(tmp)], check=True, capture_output=True)
    subprocess.run([str(SAMTOOLS), "index", str(out)], check=True)
    tmp.unlink()
    return str(out)


def first_mates(bam: str) -> str:
    r1 = Path(bam).with_suffix(".R1.bam")
    if not r1.exists():
        subprocess.run([str(SAMTOOLS), "view", "-b", "-f", "64", "-o", str(r1), bam], check=True)
        subprocess.run([str(SAMTOOLS), "index", str(r1)], check=True)
    return str(r1)


def profileseg(vals, edges, a_idx, b_idx, truth):
    smp = [dict(sample_id=f"s{i}", condition="A" if i in a_idx else "B") for i in a_idx + b_idx]
    use = [i - 1 for i in a_idx + b_idx]
    regions = [dict(chrom="c", start=k, end=k + 1, name=f"r{k}", strand="+") for k in range(vals.shape[1])]
    cfg = Config(reference="A", alternative="B", anchor="center", orientation="oriented", flank=FLANK, binsize=BIN,
                 background_width=300, core_flank=CORE)
    t0 = time.perf_counter()
    reg = per_region_tests(analyse(Profiles(vals[use], edges, smp, regions, True), cfg, bootstrap=0))
    res = {"seconds": time.perf_counter() - t0}
    for m in ("shape_kernel", "shift"):
        t = reg.tests[m]
        called = np.isin(t["call"], ("Higher", "Lower", "Changed"))
        tr = truth[t["region_index"]]
        res[m] = dict(tested=int(len(called)), called=int(called.sum()), false=int((called & ~tr).sum()),
                      true=int((called & tr).sum()), truly_changed=int(tr.sum()))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams", required=True)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--genes", type=int, default=300)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    names = [str(x) for x in np.load(args.anchors, allow_pickle=True)["anchors"]]
    anchors = [(n.split(":")[0], int(n.split(":")[1])) for n in names if n.split(":")[0] in ("chr1", "chr2")]
    rng = np.random.default_rng(13)
    subset = sorted(rng.choice(len(anchors), min(args.genes, len(anchors)), replace=False).tolist())
    anchors = [anchors[i] for i in subset]
    bed = out / "anchors.bed"
    bed.write_text("".join(f"{c}\t{max(p - FLANK - 600, 0)}\t{p + FLANK + 600}\n" for c, p in anchors))
    bams = {k: str(Path(args.bams) / f"A673_{k}.bam") for k in range(1, 7)}
    chosen_idx = np.sort(rng.choice(len(anchors), int(0.2 * len(anchors)), replace=False))
    results = {}
    script = out / "mmd.R"; script.write_text(MMD_R)
    for label, (a_idx, b_idx) in SPLITS.items():
        use_bams = dict(bams)
        truth = np.zeros(len(anchors), bool)
        if label.startswith("insertion"):
            truth[chosen_idx] = True
            for k in b_idx:
                use_bams[k] = shifted_bam(bams[k], [anchors[i] for i in chosen_idx], out / f"A673_{k}.shifted.bam")
        vals = midpoints(use_bams, anchors, bed)
        edges = np.arange(-FLANK, FLANK + BIN, BIN, dtype=float)
        res = {"ProfileSeg": profileseg(vals, edges, a_idx, b_idx, truth)}
        run_dir = out / label.split(",")[0].replace(" ", "_").replace("+", "p") / "_".join(map(str, a_idx))
        run_dir.mkdir(parents=True, exist_ok=True)
        rows = ["SampleID,Condition,Replicate,bamReads"]
        for j, k in enumerate(a_idx + b_idx):
            rows.append(f"A673_{k},{1 if k in a_idx else 2},{j % 3 + 1},{first_mates(use_bams[k])}")
        (run_dir / "samples.csv").write_text("\n".join(rows) + "\n")
        (run_dir / "peaks.tsv").write_text("chrom\tstart\tend\tname\n" + "".join(
            f"{c}\t{p - FLANK}\t{p + FLANK}\tp{i}\n" for i, (c, p) in enumerate(anchors)))
        t0 = time.perf_counter()
        proc = subprocess.run(["Rscript", str(script), str(run_dir)], capture_output=True, text=True)
        if proc.returncode:
            res["MMDiff2"] = dict(error=proc.stderr[-600:])
        else:
            mm = json.loads((run_dir / "mmd.json").read_text())
            p = np.nan_to_num(np.array(mm["p"], float), nan=1.0)
            called = bh(p) < 0.05
            res["MMDiff2"] = dict(seconds=time.perf_counter() - t0, tested=int(len(p)), called=int(called.sum()),
                                  false=int((called & ~truth).sum()), true=int((called & truth).sum()),
                                  truly_changed=int(truth.sum()), raw_p05=float(np.mean(p < 0.05)))
        results[label] = res
        print(label, json.dumps(res), flush=True)
    (out / "realbam_cutrun.json").write_text(json.dumps(results, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
