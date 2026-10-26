"""ProfileSeg versus common alternatives, on simulations with known truth.

A. Global shape change ("did the profile shape change?"); share of runs with a call:
   naive          Welch on sd of replicate aggregate profiles without background correction
   pointwise_rep  per-bin Welch t on area-normalised replicate aggregate profiles, BH over bins,
                  call if any bin q < 0.05 (the usual functional-data baseline)
   pointwise_reg  per-bin paired t with *regions* as units (condition means per region,
                  area-normalised), BH over bins: the unit choice of IWTomics / region-level curve
                  tests, which ignores replicate variability
   profileseg     any shape metric called in the report (primary + Holm)
B. Per-region amplitude, FDR and power:
   welch          ordinary per-region two-sample t on log2_signal_relative, BH
   moderated      ProfileSeg empirical-Bayes moderated t (limma-equivalent), BH
   edger_tmm      edgeR quasi-likelihood on window read counts with TMM normalisation (R), BH

Usage:  python benchmarks/profileseg_comparison.py --out comparison [--reps 40]
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import stats

from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.benchmark import region_truth
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.regional import amplitude_matrices, bh, per_region_tests, relative_offset
from genomecanvas.profileseg.simulate import SimSettings, simulate

SHAPE_SCENARIOS = ("null", "amplitude", "global_scale", "shift", "broaden", "fragment_length")
REGION_SCENARIOS = ("null", "amplitude", "global_scale", "site_amplitude", "global_with_site")
SHAPE_TRUTH = {"null": False, "amplitude": False, "global_scale": False, "shift": True, "broaden": True,
               "fragment_length": True}

EDGER = r"""
suppressMessages(library(edgeR))
args <- commandArgs(TRUE)
for (f in list.files(args[1], pattern = "^counts_.*\\.tsv$", full.names = TRUE)) {
  x <- as.matrix(read.table(f, header = FALSE))
  g <- factor(c(rep("ref", as.integer(args[2])), rep("alt", ncol(x) - as.integer(args[2]))), levels = c("ref", "alt"))
  y <- calcNormFactors(DGEList(x, group = g), method = "TMM")
  d <- model.matrix(~g)
  y <- estimateDisp(y, d)
  fit <- glmQLFit(y, d)
  res <- glmQLFTest(fit, coef = 2)
  write.table(data.frame(logFC = res$table$logFC, p = res$table$PValue), sub("counts_", "edger_", f),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
"""


def pointwise_calls(a, ref, alt, fdr=0.05):
    """Return (replicate-unit call, region-unit call) for a pointwise shape test."""
    v = a.values[:, a.complete]  # samples × regions × bins, raw (no background correction)
    agg = v.mean(axis=1)
    agg = agg / agg.sum(axis=1, keepdims=True)  # area-normalised replicate aggregates
    p_rep = stats.ttest_ind(agg[alt], agg[ref], axis=0, equal_var=False).pvalue
    by_region = [v[g].mean(axis=0) for g in (ref, alt)]
    by_region = [r / np.maximum(r.sum(axis=1, keepdims=True), M.EPS) for r in by_region]
    p_reg = stats.ttest_rel(by_region[1], by_region[0], axis=0).pvalue
    return bool(np.any(bh(p_rep) < fdr)), bool(np.any(bh(p_reg) < fdr))


def shape_comparison(reps, settings, config, seed=70_000):
    rows = []
    for k, scenario in enumerate(SHAPE_SCENARIOS):
        hits = {m: 0 for m in ("naive", "pointwise_rep", "pointwise_reg", "profileseg")}
        for rep in range(reps):
            profiles, _ = simulate(scenario, settings, seed=seed + 1000 * k + rep)
            a = analyse(profiles, config, bootstrap=50)
            n_ref = len(a.ref_idx)
            ref, alt = np.arange(n_ref), np.arange(n_ref, n_ref + len(a.alt_idx))
            eff = {e["metric"]: e for e in a.effects}
            naive = eff["naive_sd_width"]
            hits["naive"] += bool(np.isfinite(naive["welch_low"]) and (naive["welch_low"] > 0 or naive["welch_high"] < 0))
            r1, r2 = pointwise_calls(a, ref, alt)
            hits["pointwise_rep"] += r1
            hits["pointwise_reg"] += r2
            hits["profileseg"] += any(e["call"] in ("higher", "lower") for e in a.effects
                                      if e["family"] == "shape" and e["metric"] != "group_w1")
        rows.append(dict(scenario=scenario, shape_truly_changes=SHAPE_TRUTH[scenario],
                         **{m: h / reps for m, h in hits.items()}))
    return rows


def region_comparison(reps, settings, config, seed=80_000, fdr=0.05):
    have_r = shutil.which("Rscript") is not None
    work = Path(tempfile.mkdtemp(prefix="ps_cmp_"))
    runs = []
    for k, scenario in enumerate(REGION_SCENARIOS):
        for rep in range(reps):
            run_seed = seed + 1000 * k + rep
            profiles, truth = simulate(scenario, settings, seed=run_seed)
            a = analyse(profiles, config, bootstrap=10)
            res = per_region_tests(a, fdr=fdr)
            mod = res.tests["log2_signal_relative"]
            idx = mod["region_index"]
            true = region_truth(scenario, settings, config, None, run_seed)
            t_idx, t_eff = true["log2_signal_relative"]
            lookup = dict(zip(t_idx.tolist(), t_eff, strict=True))
            eff = np.array([lookup[int(i)] for i in idx])
            tr = np.where(np.abs(eff) > 0.1, np.sign(eff), 0).astype(int)
            # ordinary Welch on the same values
            n_ref = len(a.ref_idx)
            ref, alt = np.arange(n_ref), np.arange(n_ref, n_ref + len(a.alt_idx))
            _, amp = amplitude_matrices(a)
            offset, _ = relative_offset(amp["log2_signal"], ref, alt)
            rel = amp["log2_signal"] - offset
            w = stats.ttest_ind(rel[alt], rel[ref], axis=0, equal_var=False)
            welch_est = rel[alt].mean(axis=0) - rel[ref].mean(axis=0)
            runs.append(dict(scenario=scenario, rep=rep, truth=tr,
                             welch=(welch_est, bh(w.pvalue)), moderated=(mod["estimate"], mod["fdr"])))
            if have_r:
                lib = np.asarray(truth["library"])
                window = np.rint(profiles.values[:, idx, :] * lib[:, None, None]).sum(axis=-1)
                np.savetxt(work / f"counts_{k}_{rep}.tsv", window.T, fmt="%d", delimiter="\t")
    if have_r:
        script = work / "edger.R"
        script.write_text(EDGER)
        subprocess.run(["Rscript", str(script), str(work), str(settings.n_ref)], check=True, capture_output=True)
        for r in runs:
            k = REGION_SCENARIOS.index(r["scenario"])
            out = np.genfromtxt(work / f"edger_{k}_{r['rep']}.tsv", names=True, delimiter="\t")
            r["edger_tmm"] = (out["logFC"], bh(out["p"]))
    shutil.rmtree(work, ignore_errors=True)
    rows = []
    methods = ["welch", "moderated"] + (["edger_tmm"] if have_r else [])
    for scenario in REGION_SCENARIOS:
        sel = [r for r in runs if r["scenario"] == scenario]
        for m in methods:
            fdp, power, disc = [], [], []
            for r in sel:
                est, q = r[m]
                called = q < fdr
                sign = np.sign(est).astype(int)
                false = called & ((r["truth"] == 0) | (sign != r["truth"]))
                fdp.append(false.sum() / max(called.sum(), 1))
                disc.append(called.sum())
                if (r["truth"] != 0).any():
                    power.append((called & (sign == r["truth"]) & (r["truth"] != 0)).sum() / (r["truth"] != 0).sum())
            rows.append(dict(scenario=scenario, method=m, discoveries=float(np.mean(disc)), fdr=float(np.mean(fdp)),
                             power=float(np.mean(power)) if power else float("nan"),
                             true_changes=float(np.mean([(r["truth"] != 0).sum() for r in sel]))))
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="profileseg_comparison")
    parser.add_argument("--reps", type=int, default=40)
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    settings = SimSettings()
    config = Config(reference="WT", alternative="KO")
    shape = shape_comparison(args.reps, settings, config)
    region = region_comparison(args.reps, settings, config)
    (out / "comparison.json").write_text(json.dumps(dict(reps=args.reps, shape=shape, region=region), indent=1))
    lines = ["# ProfileSeg vs alternatives", "", f"{args.reps} runs per scenario, 3 vs 3, 300 regions.", "",
             "## A. Global shape change: share of runs with a call", "",
             "| scenario | shape truly changes | naive | pointwise (replicates) | pointwise (regions as units) | ProfileSeg |",
             "|---|---|---|---|---|---|"]
    for r in shape:
        lines.append(f"| {r['scenario']} | {'yes' if r['shape_truly_changes'] else 'no'} | {r['naive']:.2f} | "
                     f"{r['pointwise_rep']:.2f} | {r['pointwise_reg']:.2f} | {r['profileseg']:.2f} |")
    lines += ["", "## B. Per-region amplitude (relative / site-specific)", "",
              "| scenario | method | truly changed | discoveries | FDR | power |", "|---|---|---|---|---|---|"]
    for r in region:
        power = f"{r['power']:.2f}" if np.isfinite(r["power"]) else "—"
        lines.append(f"| {r['scenario']} | {r['method']} | {r['true_changes']:.1f} | {r['discoveries']:.1f} | "
                     f"{r['fdr']:.3f} | {power} |")
    (out / "COMPARISON.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
