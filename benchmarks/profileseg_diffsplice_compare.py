"""ProfileSeg v3 versus edgeR QL diffSpliceDGE on the same simulated counts (per-region shape changes).

diffSplice treats each region as a "gene" and its bins as "exons" and tests whether the relative
distribution of counts over bins differs between conditions (gene-level F, Simes also reported by
edgeR; the F p-value is used). It answers the same detection question as ProfileSeg's shape tests,
on counts, with a negative-binomial GLM, but has no notion of spatial order and gives no bp effect.

Per run: 600 regions, 20 % changed, the rest unchanged; counts are the simulator's integer counts
(values × library). Methods:
  ProfileSeg  shape_any, shape_kernel, shift (v3 default, full window)
  diffSplice  (edgeR QL, diffSpliceDGE) six configurations on core-window bins:
              all_F        all bins, gene-level F
              support_F    support bins (condition-blind support of the background-corrected pooled profile,
                           ≥ 20 % of its maximum, as in ProfileSeg's kernel normalisation), gene-level F
              support_simes  support bins, Simes-adjusted exon-level p
              support_agg40_F / support_agg60_F  support bins merged into 40 / 60 bp, gene-level F
              support_filter_F  support bins after filterByExpr, gene-level F
Outcomes at BH 5 % with 95 % Monte Carlo CIs, and power / FDP at BH levels 0.01–0.2 so methods can be
compared at matched observed FDP (calibration differs between methods).

Round 5: the broadening scenario's library rescaling was inverted in round 4 (fixed); support for
diffSplice now uses background-corrected pooled counts; six configurations; 50 runs per cell.

Round 6 (pre-specified before the run, seeds 9 800 000+): primary endpoint = end-to-end power at BH 5 %
(denominator: every truly changed region in the input; regions a method does not return count as not
detected); primary comparison = shape_any versus dS_support_agg60_F (the comparator chosen in round 5 and
fixed here), paired bootstrap 95 % CI of the difference over runs; hypothesis: CI > 0 for shift50_sd60 and
shift100_sd150 at 2 vs 2 and 3 vs 3. Secondary: power on the regions returned by every method (common
set), regions retained per method, full power–FDP curves over BH levels with CIs (no interpolation; the
earlier matched-FDP interpolation is dropped because mean FDP is not monotone in the BH level for
low-power curves). The run directory receives the uncommitted source diff (src.patch).
Scenarios: translations (20 / 50 / 100 / 200 bp) on narrow (sd 60) and medium (sd 150) peaks;
broadening ×1.5; a one-sided gain (half the apex, 150 bp downstream); controls without shape change:
null, global amplitude ×2, site-specific amplitude ×2 in 20 % of regions.

Outputs: power, FDP, null FWER (share of runs with any call when nothing changed in shape), and for
ProfileSeg's reliable translations the median estimate ÷ true and the direction accuracy.

    python benchmarks/profileseg_diffsplice_compare.py --out DIR [--reps 20]
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import tempfile
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.regional import bh, per_region_tests, support_mask
from genomecanvas.profileseg.simulate import SimSettings, simulate

SEED = 9_500_000
SCENARIOS = ([(f"shift{d}_sd{w}", dict(kind="shift", d=d, sigma=w)) for w in (60, 150) for d in (20, 50, 100, 200)]
             + [("broaden_x1.5", dict(kind="broaden")), ("side_gain150", dict(kind="side_gain", d=150)),
                ("null", dict(kind="null")), ("amplitude_x2", dict(kind="amplitude")), ("site_amplitude", dict(kind="site_amplitude"))])
R_SCRIPT = r"""
suppressMessages(library(edgeR))
args <- commandArgs(TRUE); run <- args[1]; n <- as.integer(args[2])
cond <- factor(c(rep("WT", n), rep("KO", n)), levels = c("WT", "KO"))
design <- model.matrix(~cond)
res <- list()
one <- function(file, filter = FALSE) {
  x <- as.matrix(read.table(file.path(run, file), header = TRUE, row.names = 1))
  region <- sub("_b[0-9]+$", "", rownames(x))
  y <- DGEList(x, genes = data.frame(GeneID = region))
  if (filter) { keep <- filterByExpr(y, design); y <- y[keep, , keep.lib.sizes = FALSE] }
  y <- normLibSizes(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  diffSpliceDGE(fit, coef = 2, geneid = "GeneID", verbose = FALSE)
}
sp <- one("counts_all.tsv"); res$all_F <- list(region = sp$gene.genes$GeneID, p = sp$gene.p.value)
sp <- one("counts_support.tsv")
res$support_F <- list(region = sp$gene.genes$GeneID, p = sp$gene.p.value)
res$support_simes <- list(region = sp$gene.genes$GeneID, p = sp$gene.Simes.p.value)
sp <- one("counts_agg40.tsv"); res$support_agg40_F <- list(region = sp$gene.genes$GeneID, p = sp$gene.p.value)
sp <- one("counts_agg60.tsv"); res$support_agg60_F <- list(region = sp$gene.genes$GeneID, p = sp$gene.p.value)
sp <- one("counts_support.tsv", filter = TRUE); res$support_filter_F <- list(region = sp$gene.genes$GeneID, p = sp$gene.p.value)
for (k in names(res)) res[[k]]$region <- as.character(res[[k]]$region)
writeLines(jsonlite::toJSON(res, digits = NA), file.path(run, "edger.json"))
"""
DSP = ("all_F", "support_F", "support_simes", "support_agg40_F", "support_agg60_F", "support_filter_F")
LEVELS = (0.01, 0.02, 0.05, 0.1, 0.2)


def make(scenario: dict, n: int, seed: int):
    kind = scenario["kind"]
    st = SimSettings(n_regions=600, n_ref=n, n_alt=n, sparse_fraction=0.2, site_fraction=0.2,
                     shift_bp=float(scenario.get("d", 100)), peak_sigma=float(scenario.get("sigma", 150)))
    rng = np.random.default_rng(seed)
    if kind == "shift":
        prof, truth = simulate("sparse_shift", st, seed=seed)
        changed = np.asarray(truth["regions"]["shift"]) != 0
    elif kind == "broaden":  # 20 % broadened: concatenate a broadened and an unchanged simulation
        a, ta = simulate("broaden", SimSettings(n_regions=120, n_ref=n, n_alt=n), seed=seed)
        b, tb = simulate("null", SimSettings(n_regions=480, n_ref=n, n_alt=n), seed=seed + 1)
        # values = counts / library: express b on a's library factors so that values × library_a returns counts_b
        # (round 4 used the inverse ratio, i.e. counts_b × (library_a / library_b)²; fixed in round 5)
        a.values = np.concatenate([a.values, b.values * (np.asarray(tb["library"]) / np.asarray(ta["library"]))[:, None, None]], axis=1)
        a.regions = a.regions + b.regions
        prof, truth = a, ta
        changed = np.r_[np.ones(120, bool), np.zeros(480, bool)]
    elif kind == "side_gain":
        prof, truth = simulate("null", st, seed=seed)
        changed = np.zeros(prof.values.shape[1], bool)
        changed[rng.choice(len(changed), int(0.2 * len(changed)), replace=False)] = True
        c = M.centers(prof.edges)
        alt = np.flatnonzero(prof.conditions == "KO")
        v = prof.values.copy()
        mean = np.nan_to_num(v).mean(0)
        lib = np.asarray(truth["library"])
        for r in np.flatnonzero(changed):
            bg = np.median(np.r_[mean[r, :10], mean[r, -10:]])
            k = int(np.argmax(mean[r]))
            bump = 0.5 * max(mean[r, k] - bg, 0) * np.exp(-0.5 * ((c - c[k] - scenario["d"]) / 40.0) ** 2)
            for s in alt:  # Poisson counts for the added signal, keeping values = counts / library
                v[s, r] += rng.poisson(bump * lib[s]) / lib[s]
        prof.values = v
    else:
        prof, truth = simulate({"null": "null", "amplitude": "amplitude", "site_amplitude": "site_amplitude"}[kind], st, seed=seed)
        changed = np.zeros(prof.values.shape[1], bool)
    lib = np.asarray(truth["library"], float)
    counts = np.rint(np.nan_to_num(prof.values) * lib[:, None, None]).astype(int)
    return prof, counts, changed


def run(job):
    name, scenario, n, rep, workdir = job[:5]
    base = job[5] if len(job) > 5 else SEED
    seed = base + 10_000 * zlib.crc32(name.encode()) % 997 * 10 + 100 * n + rep
    prof, counts, changed = make(scenario, n, seed)
    cfg = Config(reference="WT", alternative="KO")
    a = analyse(prof, cfg, bootstrap=0)
    reg = per_region_tests(a)
    out = {}

    def score(idx, called):
        t = changed[idx]
        nc = int(called.sum())
        return dict(power=float((called & t).sum() / t.sum()) if t.any() else float("nan"),
                    fdp=float((called & ~t).sum() / nc) if nc else 0.0, any_false=float((called & ~t).any()))

    returned = {}  # method → (region indices returned, p-values)
    for m in ("shape_any", "shape_kernel", "shift"):
        t = reg.tests[m]
        out[m] = score(t["region_index"], np.isin(t["call"], ("Higher", "Lower", "Changed")))
        returned[m] = (np.asarray(t["region_index"]), np.asarray(t["p"], float))
    t = reg.tests["shift"]
    if scenario["kind"] == "shift":
        idx = t["region_index"]
        rel = np.asarray(t["translation_consistent"], bool) & np.isin(t["call"], ("Higher", "Lower")) & changed[idx]
        est = np.asarray(t["estimate"], float)[rel]
        out["shift"]["reliable_share"] = float(rel.sum() / max((np.isin(t["call"], ("Higher", "Lower", "Changed")) & changed[idx]).sum(), 1))
        out["shift"]["estimate_ratio"] = float(np.median(est / scenario["d"])) if rel.any() else float("nan")
        out["shift"]["direction_ok"] = float(np.mean(est > 0)) if rel.any() else float("nan")
    # diffSplice on core-window counts of the same regions
    core = M.core_mask(prof.edges, cfg.core_flank)
    cnt = counts[..., core]
    pooled = cnt.sum(0).astype(float)
    # support from the background-corrected pooled profile (background = median of the outer tenths)
    nb = max(1, pooled.shape[1] // 10)
    bgc = np.median(np.concatenate([pooled[:, :nb], pooled[:, -nb:]], axis=1), axis=1, keepdims=True)
    sup = support_mask((pooled - bgc)[None], 0.2)
    run_dir = Path(tempfile.mkdtemp(dir=workdir))
    ids = [s["sample_id"] for s in prof.samples]
    order = np.r_[np.flatnonzero(prof.conditions == "WT"), np.flatnonzero(prof.conditions == "KO")]
    width = float(np.median(np.diff(prof.edges)))
    for kind, mask, merge in (("all", np.ones_like(sup, bool), 1), ("support", sup, 1),
                              ("agg40", sup, max(1, int(round(40 / width)))), ("agg60", sup, max(1, int(round(60 / width))))):
        rows, names = [], []
        for r in range(cnt.shape[1]):
            bins = np.flatnonzero(mask[r])
            for g, start in enumerate(range(0, len(bins), merge)):
                col = cnt[order, r][:, bins[start:start + merge]].sum(1)
                if col.sum() > 0:
                    rows.append(col); names.append(f"r{r}_b{g}")
        with open(run_dir / f"counts_{kind}.tsv", "w") as fh:
            fh.write("\t".join([ids[i] for i in order]) + "\n")
            for nm, col in zip(names, rows, strict=True):
                fh.write(nm + "\t" + "\t".join(map(str, col.tolist())) + "\n")
    subprocess.run(["Rscript", str(Path(workdir) / "edger.R"), str(run_dir), str(n)], capture_output=True, text=True, check=True)
    ed = json.loads((run_dir / "edger.json").read_text())
    for kind in DSP:
        reg_idx = np.array([int(x[1:]) for x in ed[kind]["region"]])
        p = np.array(ed[kind]["p"], float).ravel()  # Simes p-values come back as a one-column matrix
        out[f"dS_{kind}"] = score(reg_idx, bh(np.nan_to_num(p, nan=1.0)) < 0.05)
        returned[f"dS_{kind}"] = (reg_idx, p)
    # common denominators: every truly changed input region (end-to-end) and regions returned by all methods
    common = set(range(len(changed)))
    for idx, _ in returned.values():
        common &= set(idx.tolist())
    common_mask = np.zeros(len(changed), bool)
    common_mask[list(common)] = True
    for m, (idx, p) in returned.items():
        q = bh(np.nan_to_num(p, nan=1.0))
        curve = {}
        for lv in LEVELS:
            called = np.zeros(len(changed), bool)
            called[idx[q < lv]] = True
            nc = int(called.sum())
            curve[f"{lv}"] = dict(e2e_power=float((called & changed).sum() / changed.sum()) if changed.any() else float("nan"),
                                  common_power=float((called & changed & common_mask).sum() / (changed & common_mask).sum())
                                  if (changed & common_mask).any() else float("nan"),
                                  fdp=float((called & ~changed).sum() / nc) if nc else 0.0)
        out[m]["curve"] = curve
        out[m]["retained"] = int(len(idx))
        out[m]["retained_changed"] = int(changed[idx].sum())
    out["_common"] = dict(regions=int(common_mask.sum()), changed=int((changed & common_mask).sum()), total_changed=int(changed.sum()))
    for f in run_dir.iterdir():
        f.unlink()
    run_dir.rmdir()
    return name, n, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--seed-base", type=int, default=SEED)
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    git = lambda *a: subprocess.run(["git", *a], cwd=root, capture_output=True, text=True).stdout.strip()  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                regional_py_sha256=hashlib.sha256((root / "src/genomecanvas/profileseg/regional.py").read_bytes()).hexdigest()[:16],
                src_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--", "src").encode()).hexdigest()[:16],
                edgeR=subprocess.run(["Rscript", "-e", 'cat(as.character(packageVersion("edgeR")))'], capture_output=True, text=True).stdout.strip(),
                reps=args.reps, seed_base=args.seed_base, started=datetime.datetime.now().isoformat(timespec="seconds"))
    (out / "src.patch").write_text(git("diff", "HEAD", "--", "src") + "\n")
    work = Path(tempfile.mkdtemp(prefix="dsp_", dir=out))
    (work / "edger.R").write_text(R_SCRIPT)
    jobs = [(nm, sc, n, r, str(work), args.seed_base) for nm, sc in SCENARIOS for n in (2, 3) for r in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(args.jobs) as ex:
        for nm, n, res in ex.map(run, jobs):
            agg.setdefault((nm, n), []).append(res)
    (work / "edger.R").unlink()
    work.rmdir()
    prov["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    methods = ("shape_any", "shape_kernel", "shift") + tuple(f"dS_{k}" for k in DSP)
    rng = np.random.default_rng(1)

    def mci(v):
        v = np.asarray(v, float); v = v[np.isfinite(v)]
        return (float(v.mean()), float(1.96 * v.std(ddof=1) / np.sqrt(len(v)))) if len(v) > 1 else (float("nan"), float("nan"))

    def paired(a, b, B=4000):
        a, b = np.asarray(a, float), np.asarray(b, float)
        d = a - b
        bs = d[rng.integers(0, len(d), (B, len(d)))].mean(1)
        return float(d.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))

    cells = {}
    for (nm, n), rows in agg.items():
        c = {}
        for m in methods:
            d = {}
            for k in ("any_false", "reliable_share", "estimate_ratio", "direction_ok", "retained", "retained_changed"):
                v = [r[m].get(k, np.nan) for r in rows]
                d[k], d[k + "_ci"] = mci(v)
            d["curve"] = {}
            for lv in rows[0][m]["curve"]:
                d["curve"][lv] = {}
                for k in ("e2e_power", "common_power", "fdp"):
                    mu, ci = mci([r[m]["curve"][lv][k] for r in rows])
                    d["curve"][lv][k], d["curve"][lv][k + "_ci"] = mu, ci
            c[m] = d
        c["_common"] = {k: float(np.mean([r["_common"][k] for r in rows])) for k in rows[0]["_common"]}
        if nm not in ("null", "amplitude_x2", "site_amplitude"):
            c["_primary"] = dict(zip(("diff", "low", "high"), paired([r["shape_any"]["curve"]["0.05"]["e2e_power"] for r in rows],
                                                                     [r["dS_support_agg60_F"]["curve"]["0.05"]["e2e_power"] for r in rows])))
        cells[f"{nm}|{n}v{n}"] = c
    (out / "diffsplice_compare.json").write_text(json.dumps(dict(provenance=prov, cells=cells), indent=1))
    hyp = [("shift50_sd60", 2), ("shift50_sd60", 3), ("shift100_sd150", 2), ("shift100_sd150", 3)]
    lines = ["# ProfileSeg v3 vs edgeR QL diffSpliceDGE (same simulated counts)", "",
             f"Commit `{prov['commit'][:10]}`; uncommitted src files {prov['uncommitted_src']}; regional.py `{prov['regional_py_sha256']}`; "
             f"src diff `{prov['src_diff_sha256']}` (patch: src.patch); edgeR {prov['edgeR']}; {args.reps} runs per cell; seeds {args.seed_base}+; "
             f"{prov['started']} → {prov['finished']}.", "",
             "## Primary (pre-specified): end-to-end power at BH 5 %, shape_any − dS_support_agg60_F", "",
             "| scenario | design | shape_any | dS_support_agg60_F | difference [paired bootstrap 95 % CI] | hypothesis |", "|---|---|---|---|---|---|"]
    for nm, _ in SCENARIOS:
        if nm in ("null", "amplitude_x2", "site_amplitude"):
            continue
        for n in (2, 3):
            c = cells[f"{nm}|{n}v{n}"]
            a5, d5 = c["shape_any"]["curve"]["0.05"], c["dS_support_agg60_F"]["curve"]["0.05"]
            pr = c["_primary"]
            h = ("**met**" if pr["low"] > 0 else "**not met**") if (nm, n) in hyp else ""
            lines.append(f"| {nm} | {n}v{n} | {a5['e2e_power']:.3f}±{a5['e2e_power_ci']:.3f} | {d5['e2e_power']:.3f}±{d5['e2e_power_ci']:.3f} | "
                         f"{pr['diff']:+.3f} [{pr['low']:+.3f}, {pr['high']:+.3f}] | {h} |")
    lines += ["", "## All methods at BH 5 %: end-to-end power / common-set power / FDP (± 95 % CI); regions retained (changed retained)", "",
              "| scenario | design | " + " | ".join(methods) + " |", "|---|---|" + "---|" * len(methods)]
    for nm, _ in SCENARIOS:
        for n in (2, 3):
            c = cells[f"{nm}|{n}v{n}"]
            if nm in ("null", "amplitude_x2", "site_amplitude"):
                vals = [f"any call {c[m]['any_false']:.2f}" for m in methods]
            else:
                vals = []
                for m in methods:
                    q = c[m]["curve"]["0.05"]
                    vals.append(f"{q['e2e_power']:.2f} / {q['common_power']:.2f} / {q['fdp']:.3f}±{q['fdp_ci']:.3f} "
                                f"({c[m]['retained']:.0f}; {c[m]['retained_changed']:.0f})")
            lines.append(f"| {nm} | {n}v{n} | " + " | ".join(vals) + " |")
    lines += ["", f"Common set: regions returned by every method (mean over runs); see JSON `_common`.", "",
              "## Power–FDP curves (end-to-end power / FDP at BH levels " + ", ".join(str(x) for x in LEVELS) + ")", "",
              "| scenario | design | method | " + " | ".join(f"BH {x}" for x in LEVELS) + " |", "|---|---|---|" + "---|" * len(LEVELS)]
    for nm, _ in SCENARIOS:
        if nm in ("null", "amplitude_x2", "site_amplitude"):
            continue
        for n in (2, 3):
            c = cells[f"{nm}|{n}v{n}"]
            for m in ("shape_any", "shape_kernel", "dS_support_agg60_F", "dS_support_F"):
                pts = [c[m]["curve"][f"{lv}"] for lv in LEVELS]
                lines.append(f"| {nm} | {n}v{n} | {m} | " + " | ".join(f"{q['e2e_power']:.2f}±{q['e2e_power_ci']:.2f} / {q['fdp']:.3f}" for q in pts) + " |")
    (out / "DIFFSPLICE_COMPARE.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
