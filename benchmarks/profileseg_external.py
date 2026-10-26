"""ProfileSeg versus IWTomics and MMDiff2 (Bioconductor) on simulations with known truth.

Both tools are run as published, through R:

* IWTomics (Cremona et al. 2018; interval-wise testing, Pini & Vantini 2017): two-sample test of two
  groups of curves, permutation-based, p-values adjusted over all intervals of the domain. Run in the
  two ways a user can apply it to a WT-vs-KO design:
    - ``iwt_regions``: curves = regions (condition-mean profile per region), paired across conditions;
    - ``iwt_replicates``: curves = replicate aggregate profiles (3 vs 3; only 20 distinct permutations).
  Each in raw units and with every curve scaled to unit area (shape only). Call = any position with
  adjusted p < 0.05.
* MMDiff2 (Schweikert et al. 2013): per-peak kernel MMD between read-centre distributions, null from
  within-condition replicate distances (locfit over coverage), on BAM files written from the simulated
  counts. Per-region calls at BH FDR < 0.05 on its p-values.

ProfileSeg answers the same two questions with its aggregate shape metrics (global) and per-region
centroid / mad_width tests (per region).

Usage:  python benchmarks/profileseg_external.py --out external [--reps 10]
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.regional import bh, per_region_tests
from genomecanvas.profileseg.simulate import SimSettings, simulate

SCENARIOS = ("null", "amplitude", "global_scale", "shift", "broaden", "opposite_shift", "site_amplitude")
SHAPE_CHANGES = {"null": False, "amplitude": False, "global_scale": False, "shift": True, "broaden": True,
                 "opposite_shift": True, "site_amplitude": False}
CORE = 1000  # bp each side; IWTomics curves and MMDiff2 peaks use the same core window as ProfileSeg
IWT_BIN = 40  # bp; coarser grid keeps interval-wise testing tractable
FRAGMENT = 200
READ = 36

R_SCRIPT = r"""
suppressMessages({library(IWTomics); library(MMDiff2); library(GenomicRanges); library(Rsamtools)})
run <- commandArgs(TRUE)[1]
iwt_call <- function(mats, paired) {
  # mats: list(group1 = bins x curves, group2 = bins x curves)
  nb <- nrow(mats[[1]])
  mk <- function(n, id) GRanges("chrSim", IRanges(seq_len(n) * 10000, width = nb * 40))
  regs <- GRangesList(g1 = mk(ncol(mats[[1]])), g2 = mk(ncol(mats[[2]])))
  feats <- list(f = list(g1 = mats[[1]], g2 = mats[[2]]))
  d <- IWTomicsData(regs, feats, alignment = "center", id_regions = c("g1", "g2"),
                    id_features = "f", length_features = list(f = list(g1 = rep(nb, ncol(mats[[1]])),
                                                                         g2 = rep(nb, ncol(mats[[2]])))))
  r <- suppressMessages(IWTomicsTest(d, id_region1 = "g2", id_region2 = "g1", paired = paired, B = 1000))
  p <- adjusted_pval(r)[[1]][[1]]
  min(p, na.rm = TRUE)
}
{
  res <- list()
  for (kind in c("raw", "area")) {
    reg <- lapply(c("ref", "alt"), function(g) as.matrix(read.table(file.path(run, sprintf("iwt_regions_%s_%s.tsv", kind, g)))))
    rep <- lapply(c("ref", "alt"), function(g) as.matrix(read.table(file.path(run, sprintf("iwt_reps_%s_%s.tsv", kind, g)))))
    res[[paste0("iwt_regions_", kind)]] <- iwt_call(reg, TRUE)
    res[[paste0("iwt_replicates_", kind)]] <- iwt_call(rep, FALSE)
  }
  # MMDiff2
  for (s in list.files(run, pattern = "\\.sam$", full.names = TRUE))
    asBam(s, sub("\\.sam$", "", s), overwrite = TRUE, indexDestination = TRUE)
  meta <- list(ExpData = list(dataDir = run, sampleSheet = "samples.csv"))
  peaks <- read.table(file.path(run, "peaks.tsv"), header = TRUE)
  gr <- GRanges(peaks$chrom, IRanges(peaks$start + 1, peaks$end), names = peaks$name)
  names(gr) <- peaks$name
  md <- DBAmmd(meta)
  md <- setRegions(md, gr)
  md <- suppressMessages(getPeakReads(md, pairedEnd = FALSE, run.parallel = FALSE))
  md <- suppressMessages(estimateFragmentCenters(md))
  md <- compHists(md, bin.length = 20, whichPos = "Center")
  md <- suppressMessages(compDists(md, dist.method = "MMD", run.parallel = FALSE))
  md <- setContrast(md, contrast = "byCondition")
  md <- suppressMessages(compPvals(md, dist.method = "MMD"))
  tab <- md@Contrasts[[1]]$MMD.locfit$mmd.Table
  res$mmdiff_p <- as.numeric(tab[, "pval"])
  res$mmdiff_padj <- as.numeric(tab[, "padj"])  # MMDiff2's own adjustment (reportResults default th = 0.1)
  res$mmdiff_outlier <- as.logical(md@Contrasts[[1]]$MMD.locfit$outlier)
  for (f in list.files(run, pattern = "\\.(sam|bam|bai)$", full.names = TRUE)) file.remove(f)
}
writeLines(jsonlite::toJSON(res, digits = NA, auto_unbox = TRUE), file.path(run, "r_result.json"))
"""


def _iwt_matrices(a, ref, alt):
    """Core-window curves on a coarser grid: per-region condition means and replicate aggregates."""
    e = a.profiles.edges
    c = M.centers(e)
    keep = np.abs(c) < CORE
    v = a.values[:, a.complete][:, :, keep]  # samples × regions × bins (raw coverage)
    k = IWT_BIN // int(round(np.diff(e)[0]))
    v = v[..., : (v.shape[-1] // k) * k].reshape(*v.shape[:2], -1, k).mean(axis=-1)
    area = lambda x: x / np.maximum(x.sum(axis=-1, keepdims=True), M.EPS)  # noqa: E731
    regions = {g: v[idx].mean(axis=0) for g, idx in (("ref", ref), ("alt", alt))}  # regions × bins
    reps = {g: v[idx].mean(axis=1) for g, idx in (("ref", ref), ("alt", alt))}  # replicates × bins
    return {
        "regions_raw": regions, "regions_area": {g: area(x) for g, x in regions.items()},
        "reps_raw": reps, "reps_area": {g: area(x) for g, x in reps.items()},
    }


def _write_sam(path: Path, counts: np.ndarray, starts: np.ndarray, edges: np.ndarray, rng, chrom_len: int):
    """Single-end reads whose fragment centres fall uniformly inside each count bin."""
    lines = [f"@HD\tVN:1.6\tSO:unsorted", f"@SQ\tSN:chrSim\tLN:{chrom_len}"]
    width = np.diff(edges)
    n = 0
    for r, anchor in enumerate(starts):
        for b in np.flatnonzero(counts[r] > 0):
            k = int(counts[r, b])
            centre = anchor + edges[b] + rng.uniform(0, width[b], k)
            plus = rng.random(k) < 0.5
            pos = np.where(plus, centre - FRAGMENT / 2, centre + FRAGMENT / 2 - READ).astype(int) + 1
            for p, s in zip(pos, plus, strict=True):
                if p < 1:
                    continue
                n += 1
                lines.append(f"r{n}\t{0 if s else 16}\tchrSim\t{p}\t60\t{READ}M\t*\t0\t0\t*\t*")
    path.write_text("\n".join(lines) + "\n")


def prepare_run(run_dir: Path, scenario: str, settings: SimSettings, config: Config, seed: int) -> dict:
    run_dir.mkdir(parents=True)
    profiles, truth = simulate(scenario, settings, seed=seed)
    a = analyse(profiles, config, bootstrap=20)
    n_ref = len(a.ref_idx)
    ref, alt = np.arange(n_ref), np.arange(n_ref, n_ref + len(a.alt_idx))
    for name, groups in _iwt_matrices(a, ref, alt).items():
        for g, x in groups.items():
            np.savetxt(run_dir / f"iwt_{name}_{g}.tsv", x.T, delimiter="\t")  # bins × curves
    rng = np.random.default_rng(seed)
    lib = np.asarray(truth["library"])
    anchors = np.array([r["start"] for r in profiles.regions])
    rows = ["SampleID,Condition,Replicate,bamReads"]
    for i, s in enumerate(profiles.samples):
        counts = np.rint(np.nan_to_num(profiles.values[i]) * lib[i]).astype(int)
        _write_sam(run_dir / f"{s['sample_id']}.sam", counts, anchors, profiles.edges, rng,
                   chrom_len=int(anchors.max() + 20_000))
        rows.append(f"{s['sample_id']},{1 if s['condition'] == 'WT' else 2},{s['biological_replicate']},{s['sample_id']}.bam")
    (run_dir / "samples.csv").write_text("\n".join(rows) + "\n")
    comp = np.flatnonzero(a.complete)
    peaks = ["chrom\tstart\tend\tname"] + [f"chrSim\t{anchors[k] - CORE}\t{anchors[k] + CORE}\tp{k}" for k in comp]
    (run_dir / "peaks.tsv").write_text("\n".join(peaks) + "\n")
    # ProfileSeg on the same data
    shape_called = any(e["call"] in ("higher", "lower") for e in a.effects
                       if e["family"] == "shape" and e["metric"] != "group_w1")
    reg = per_region_tests(a)
    per_region = {}
    for m in ("centroid", "mad_width"):
        t = reg.tests[m]
        per_region[m] = dict(index=t["region_index"].tolist(), called=np.isin(t["call"], ("Higher", "Lower")).tolist())
    return dict(scenario=scenario, profileseg_global=bool(shape_called), profileseg_regions=per_region,
                complete=comp.tolist())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="profileseg_external")
    parser.add_argument("--reps", type=int, default=10)
    parser.add_argument("--regions", type=int, default=200)
    parser.add_argument("--jobs", type=int, default=8, help="parallel R processes")
    parser.add_argument("--keep", action="store_true", help="keep the BAM/matrix working directory")
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    settings = SimSettings(n_regions=args.regions)
    config = Config(reference="WT", alternative="KO")
    work = Path(tempfile.mkdtemp(prefix="ps_ext_", dir=out))
    runs = {}
    for k, scenario in enumerate(SCENARIOS):
        for rep in range(args.reps):
            key = f"{scenario}__{rep}"
            runs[key] = prepare_run(work / key, scenario, settings, config, seed=60_000 + 1000 * k + rep)
    script = work / "external.R"
    script.write_text(R_SCRIPT)

    def run_r(key):
        proc = subprocess.run(["Rscript", str(script), str(work / key)], capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError(f"{key}: {proc.stderr[-2000:]}")
        print(f"  R done: {key}", flush=True)
        return key, json.loads((work / key / "r_result.json").read_text())

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        r = dict(pool.map(run_r, runs))
    methods = ["iwt_regions_raw", "iwt_regions_area", "iwt_replicates_raw", "iwt_replicates_area"]
    table, regional = [], []
    # MMDiff2's own p-values are anti-conservative (Fisher combination of dependent pairwise p-values);
    # "calibrated" uses the raw-p threshold that calls 5 % of regions in the null runs.
    null_p = np.concatenate([np.asarray(r[k]["mmdiff_p"], dtype=float) for k in runs if runs[k]["scenario"] == "null"])
    threshold = float(np.nanquantile(null_p, 0.05))
    for scenario in SCENARIOS:
        keys = [k for k in runs if runs[k]["scenario"] == scenario]
        row = dict(scenario=scenario, shape_changes=SHAPE_CHANGES[scenario],
                   profileseg=float(np.mean([runs[k]["profileseg_global"] for k in keys])))
        for m in methods:
            row[m] = float(np.mean([r[k][m] < 0.05 for k in keys]))
        row["iwt_replicates_min_p"] = float(np.min([r[k]["iwt_replicates_area"] for k in keys]))
        table.append(row)
        # per region: share of regions called (truth: every region changes shape in shift/broaden/opposite)
        frac = {"mmdiff2": [], "mmdiff2_own": [], "mmdiff2_calibrated": [], "profileseg_centroid": [],
                "profileseg_mad_width": [], "profileseg_any": []}
        for k in keys:
            p = np.asarray(r[k]["mmdiff_p"], dtype=float)
            frac["mmdiff2"].append(float(np.mean(np.nan_to_num(bh(p), nan=1) < 0.05)))
            own = np.asarray(r[k]["mmdiff_padj"], dtype=float)
            frac["mmdiff2_own"].append(float(np.mean(np.nan_to_num(own, nan=1) < 0.1)))
            frac["mmdiff2_calibrated"].append(float(np.mean(np.nan_to_num(p, nan=1) <= threshold)))
            pr = runs[k]["profileseg_regions"]
            n_regions = len(runs[k]["complete"])
            called = set()
            for m in ("centroid", "mad_width"):
                c = {i for i, f in zip(pr[m]["index"], pr[m]["called"], strict=True) if f}
                frac[f"profileseg_{m}"].append(len(c) / n_regions)
                called |= c
            frac["profileseg_any"].append(len(called) / n_regions)
        regional.append(dict(scenario=scenario, shape_changes=SHAPE_CHANGES[scenario],
                             **{m: float(np.nanmean(v)) for m, v in frac.items()}))
    (out / "external.json").write_text(json.dumps(dict(reps=args.reps, regions=args.regions, global_=table,
                                                       regional=regional), indent=1))
    lines = ["# ProfileSeg vs IWTomics and MMDiff2", "",
             f"{args.reps} runs per scenario, 3 vs 3, {args.regions} regions, core ±{CORE} bp.", "",
             "## Global: share of runs reporting a profile/shape difference", "",
             "| scenario | shape changes | IWTomics regions, raw | IWTomics regions, area | IWTomics replicates, raw | "
             "IWTomics replicates, area | ProfileSeg |", "|---|---|---|---|---|---|---|"]
    for t in table:
        lines.append(f"| {t['scenario']} | {'yes' if t['shape_changes'] else 'no'} | {t['iwt_regions_raw']:.2f} | "
                     f"{t['iwt_regions_area']:.2f} | {t['iwt_replicates_raw']:.2f} | {t['iwt_replicates_area']:.2f} | "
                     f"{t['profileseg']:.2f} |")
    lines += ["", f"Smallest IWTomics adjusted p with replicates as curves: "
              f"{min(t['iwt_replicates_min_p'] for t in table):.3f}.", "",
              "## Per region: share of regions called (BH FDR < 0.05)", "",
              f"MMDiff2 calibrated: raw p ≤ {threshold:.3g}, the threshold that calls 5 % of regions in the null runs.", "",
              "| scenario | shape changes | MMDiff2 (BH < 0.05) | MMDiff2 (own padj < 0.1) | MMDiff2 calibrated | "
              "ProfileSeg centroid | ProfileSeg mad_width | ProfileSeg either |",
              "|---|---|---|---|---|---|---|---|"]
    for t in regional:
        lines.append(f"| {t['scenario']} | {'yes (all regions)' if t['shape_changes'] else 'no'} | {t['mmdiff2']:.3f} | "
                     f"{t['mmdiff2_own']:.3f} | {t['mmdiff2_calibrated']:.3f} | "
                     f"{t['profileseg_centroid']:.3f} | {t['profileseg_mad_width']:.3f} | {t['profileseg_any']:.3f} |")
    (out / "EXTERNAL.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
