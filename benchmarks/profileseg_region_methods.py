"""Per-region shape inference: candidate methods compared on the same simulated data.

Candidates (all keep biological replicates as the unit of inference):

  A0  current       centroid / MAD over the full core window, moderated t
  A   support       centroid / MAD restricted to the condition-blind support of the pooled profile
  B   weighted      A with per-sample precision weights (delta-method variance from coverage)
  C1  count GLM     edgeR QL on bin counts, diffSpliceDGE gene-level F (region = gene, bin = exon)
  C2  template fit  per-sample shift / log-scale from a weighted Gauss–Newton fit to the pooled
                    template (efficient, unbiased effect sizes), moderated t
  D   group test    per-sample mean of C2 over groups of 20 regions, moderated t over groups
  E   kernel test   kernel (MMD-type) statistic on normalised support profiles, null pooled over
                    regions and replicate-label permutations

Scenarios: the simulator's whole-set scenarios plus 20 % mixtures (60 changed regions among 300),
which is the realistic case for per-region testing. Usage:

    python benchmarks/profileseg_region_methods.py --out DIR [--templates real_templates.npz] [--reps 10]
"""

from __future__ import annotations

import argparse
import os
import itertools
import json
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

# One BLAS thread per worker process; nested BLAS threads stalled the thread-pool version.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.regional import bh, fit_f_dist, moderated_t
from genomecanvas.profileseg.simulate import SimSettings, simulate

FDR = 0.05
MIX_FRACTION = 0.2
GROUP = 20
COUNT_BIN = 3  # 20 bp simulation bins merged to 60 bp for the count GLM
SCENARIOS = ("null", "amplitude", "site_amplitude", "shift", "broaden", "mix_shift", "mix_broaden")
TRUTH = {"null": None, "amplitude": None, "site_amplitude": None, "shift": "shift", "broaden": "width",
         "mix_shift": "shift", "mix_broaden": "width"}

R_SCRIPT = r"""
suppressMessages(library(edgeR))
run <- commandArgs(TRUE)[1]
res <- list()
for (kind in c("all", "support")) {
  f <- file.path(run, sprintf("counts_%s.tsv", kind))
  x <- as.matrix(read.table(f, header = TRUE, row.names = 1))
  region <- sub("_b[0-9]+$", "", rownames(x))
  cond <- factor(c(rep("WT", 3), rep("KO", 3)), levels = c("WT", "KO"))
  design <- model.matrix(~cond)
  y <- DGEList(x)
  y <- normLibSizes(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  sp <- diffSpliceDGE(fit, coef = 2, geneid = region, verbose = FALSE)
  res[[kind]] <- list(region = as.character(sp$gene.genes$GeneID), p = as.numeric(sp$gene.p.value))
}
writeLines(jsonlite::toJSON(res, digits = NA), file.path(run, "edger.json"))
"""


@dataclass
class Data:
    y: np.ndarray        # samples × regions × bins, background-corrected core coverage
    cov: np.ndarray      # samples × regions × bins, raw core coverage (y + background)
    counts: np.ndarray   # samples × regions × bins, integer counts in the core
    centers: np.ndarray  # bin centres (bp)
    changed: np.ndarray  # regions × bool (truth for this run)


def prepare(scenario: str, settings: SimSettings, templates, seed: int) -> Data:
    config = Config(reference="WT", alternative="KO")

    def one(name, s):
        prof, truth = simulate(name, settings, seed=s, templates=templates)
        a = analyse(prof, config, bootstrap=0)
        ok = a.complete & (np.nanmean(a.core_mass, axis=0) > 0)
        core = M.core_mask(prof.edges, config.core_flank)
        lib = np.asarray(truth["library"])[:, None, None]
        cov = np.nan_to_num(prof.values[:, :, core])
        counts = np.rint(cov * lib).astype(int)
        return Data(np.nan_to_num(a.core_signal)[:, ok], cov[:, ok], counts[:, ok], M.centers(a.core_edges),
                    np.full(int(ok.sum()), TRUTH[name] is not None))

    if not scenario.startswith("mix_"):
        return one(scenario, seed)
    changed, base = one(scenario[4:], seed), one("null", seed + 50_000)
    n = min(changed.y.shape[1], base.y.shape[1])
    k = int(round(MIX_FRACTION * n))
    cat = lambda f: np.concatenate([getattr(changed, f)[:, :k], getattr(base, f)[:, : n - k]], axis=1)
    return Data(cat("y"), cat("cov"), cat("counts"), changed.centers,
                np.r_[np.ones(k, bool), np.zeros(n - k, bool)])


REF, ALT = np.arange(3), np.arange(3, 6)


def support(d: Data, frac: float = 0.2) -> np.ndarray:
    pooled = d.y.mean(axis=0)
    return pooled > frac * pooled.max(axis=1, keepdims=True)


def moments(d: Data, mask: np.ndarray | None):
    y = d.y if mask is None else np.clip(d.y, 0, None) * mask
    w = np.ones_like(d.centers)
    p = M.normalise(y * w)
    edges = np.r_[d.centers - 10, d.centers[-1] + 10]
    return M.centroid(p, edges), M.mean_abs_deviation(p, edges)


def weighted_moderated_t(v: np.ndarray, var: np.ndarray) -> dict:
    """limma-style weighted two-group fit: weights 1/var (normalised within region), eBayes on s²."""
    w = 1.0 / np.maximum(var, 1e-12)
    w = w / w.mean(axis=0, keepdims=True)
    w1, w2 = w[REF], w[ALT]
    m1 = (w1 * v[REF]).sum(0) / w1.sum(0)
    m2 = (w2 * v[ALT]).sum(0) / w2.sum(0)
    s2 = ((w1 * (v[REF] - m1) ** 2).sum(0) + (w2 * (v[ALT] - m2) ** 2).sum(0)) / 4
    d0, s02 = fit_f_dist(s2, 4)
    post = s2 if not np.isfinite(d0) else (d0 * s02 + 4 * s2) / (d0 + 4)
    t = (m2 - m1) / np.sqrt(post * (1 / w1.sum(0) + 1 / w2.sum(0)))
    from scipy import stats
    df = 4 + (d0 if np.isfinite(d0) else 1e6)
    return dict(estimate=m2 - m1, p=2 * stats.t.sf(np.abs(t), df))


def centroid_variance(d: Data, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Delta-method variance of the support centroid and MAD; bin variance ∝ raw coverage (Poisson)."""
    y = np.clip(d.y, 0, None) * mask
    mass = y.sum(-1, keepdims=True) + 1e-12
    mu = (y * d.centers).sum(-1, keepdims=True) / mass
    var_bin = d.cov * mask
    vc = (((d.centers - mu) ** 2) * var_bin).sum(-1) / mass[..., 0] ** 2
    mad = (y * np.abs(d.centers - mu)).sum(-1, keepdims=True) / mass
    vm = (((np.abs(d.centers - mu) - mad) ** 2) * var_bin).sum(-1) / mass[..., 0] ** 2
    return vc, vm


def template_fit(d: Data, iterations: int = 4, smooth: int = 5):
    """Per-sample shift (bp) and log2 scale relative to the pooled template, by weighted Gauss–Newton.

    Model per sample: y(x) ≈ k·T((x − μ − δ)/g + μ). Linearised around the current fit, the
    columns T, T' and (x − μ)·T' give the amplitude, shift and dilation updates. Weights are
    1/coverage (Poisson). Bins where the template is below 5 % of its maximum are left out.
    """
    c = d.centers
    T0 = np.clip(d.y.mean(axis=0), 0, None)
    if smooth > 1:  # light triangular smoothing of the pooled template
        h = (smooth + 1) // 2
        kernel = np.r_[np.arange(1, h + 1), np.arange(h - 1, 0, -1)].astype(float)
        kernel /= kernel.sum()
        T0 = np.apply_along_axis(lambda r: np.convolve(r, kernel, mode="same"), 1, T0)
    mu = (T0 * c).sum(1) / (T0.sum(1) + 1e-12)
    use = T0 > 0.05 * T0.max(1, keepdims=True)
    S, R, _ = d.y.shape
    delta, logg = np.zeros((S, R)), np.zeros((S, R))
    var_d, var_g = np.full((S, R), np.nan), np.full((S, R), np.nan)
    for r in range(R):
        u = use[r]
        for s in range(S):
            dl, lg = 0.0, 0.0
            for _ in range(iterations):
                g = 2.0 ** lg
                xs = (c - mu[r] - dl) / g + mu[r]
                T = np.interp(xs, c, T0[r], left=0, right=0)
                dT = np.gradient(T, c)
                X = np.column_stack([T, -dT / g, -(xs - mu[r]) * dT * np.log(2)])[u]
                wt = 1.0 / np.maximum(d.cov[s, r, u], 1.0)
                yy = d.y[s, r, u]
                XtW = X.T * wt
                A = XtW @ X
                try:
                    beta = np.linalg.solve(A, XtW @ yy)
                except np.linalg.LinAlgError:
                    break
                k = beta[0] if abs(beta[0]) > 1e-9 else 1e-9
                dl = float(np.clip(dl + beta[1] / k, -400, 400))
                lg = float(np.clip(lg + beta[2] / k, -1.5, 1.5))
            delta[s, r], logg[s, r] = dl, lg
            try:
                resid = yy - X @ beta
                sigma2 = float((wt * resid**2).sum() / max(u.sum() - 3, 1))
                cov = np.linalg.inv(A) * sigma2
                var_d[s, r] = cov[1, 1] / k**2
                var_g[s, r] = cov[2, 2] / k**2
            except (np.linalg.LinAlgError, UnboundLocalError):
                pass
    return delta, logg, var_d, var_g


def kernel_test(d: Data, mask: np.ndarray, sigma_bp: float = 60.0, s0_quantile: float = 0.5) -> np.ndarray:
    y = np.clip(d.y, 0, None) * mask
    p = y / (y.sum(-1, keepdims=True) + 1e-12)
    c = d.centers
    K = np.exp(-0.5 * ((c[:, None] - c[None, :]) / sigma_bp) ** 2)
    splits = [s for s in itertools.combinations(range(6), 3) if 0 in s]  # 10 distinct splits
    F = []
    for s in splits:
        a = np.array(s); b = np.setdiff1d(np.arange(6), a)
        ma, mb = p[a].mean(0), p[b].mean(0)
        diff = ma - mb
        D = np.einsum("rb,bc,rc->r", diff, K, diff)
        res = np.concatenate([p[a] - ma, p[b] - mb])
        W = np.einsum("srb,bc,src->r", res, K, res) / 4
        F.append(D / (W + np.quantile(W, s0_quantile)))
    F = np.array(F)  # split 0 = (0,1,2) vs (3,4,5) = observed
    obs, null = F[0], F[1:].ravel()
    null = np.sort(null)
    return (1 + null.size - np.searchsorted(null, obs, side="left")) / (1 + null.size)


def count_tables(d: Data, mask: np.ndarray, run_dir: Path) -> None:
    S, R, B = d.counts.shape
    nb = B // COUNT_BIN
    agg = d.counts[..., : nb * COUNT_BIN].reshape(S, R, nb, COUNT_BIN).sum(-1)
    m = mask[:, : nb * COUNT_BIN].reshape(R, nb, COUNT_BIN).any(-1)
    for kind, keep in (("all", np.ones_like(m)), ("support", m)):
        rows = ["\t".join(["id"] + [f"s{i}" for i in range(S)])]
        for r in range(R):
            for b in np.flatnonzero(keep[r]):
                rows.append("\t".join([f"r{r}_b{b}"] + [str(int(v)) for v in agg[:, r, b]]))
        (run_dir / f"counts_{kind}.tsv").write_text("\n".join(rows) + "\n")


def evaluate(d: Data, run_dir: Path) -> dict:
    mask = support(d)
    out = {}
    cf, mf = moments(d, None)
    out["A0_centroid"] = moderated_t(cf, REF, ALT)
    out["A0_mad"] = moderated_t(mf, REF, ALT)
    cs, ms = moments(d, mask)
    out["A_centroid"] = moderated_t(cs, REF, ALT)
    out["A_mad"] = moderated_t(ms, REF, ALT)
    vc, vm = centroid_variance(d, mask)
    out["B_centroid"] = weighted_moderated_t(cs, vc)
    out["B_mad"] = weighted_moderated_t(ms, vm)
    delta, logg, vd, vg = template_fit(d)
    out["C2_shift"] = moderated_t(delta, REF, ALT)
    out["C2_scale"] = moderated_t(logg, REF, ALT)
    out["C2w_shift"] = weighted_moderated_t(delta, np.nan_to_num(vd, nan=np.nanmax(vd)))
    out["C2w_scale"] = weighted_moderated_t(logg, np.nan_to_num(vg, nan=np.nanmax(vg)))
    out["E_kernel"] = dict(p=kernel_test(d, mask), estimate=np.full(d.y.shape[1], np.nan))
    count_tables(d, mask, run_dir)
    # D: groups of GROUP consecutive regions (mixtures put the changed regions first)
    R = d.y.shape[1]
    ng = R // GROUP
    for name, v in (("D_shift", delta), ("D_scale", logg)):
        gm = v[:, : ng * GROUP].reshape(6, ng, GROUP).mean(-1)
        out[name] = moderated_t(gm, REF, ALT)
    return {k: dict(p=np.asarray(v["p"], float).tolist(), estimate=np.asarray(v["estimate"], float).tolist())
            for k, v in out.items()}, d.changed.tolist(), ng


def score(p, truth) -> dict:
    p = np.nan_to_num(np.asarray(p, float), nan=1.0)
    called = bh(p) < FDR
    truth = np.asarray(truth, bool)
    n_called = int(called.sum())
    return dict(called=float(called.mean()), any=float(called.any()),
                fdp=float((called & ~truth).sum() / n_called) if n_called else 0.0,
                power=float((called & truth).sum() / truth.sum()) if truth.any() else float("nan"))


def run_one(job):
    sc, rep, work, settings, templates = job
    run_dir = work / f"{sc}__{rep}"; run_dir.mkdir()
    d = prepare(sc, settings, templates, seed=70_000 + 1000 * SCENARIOS.index(sc) + rep)
    res, truth, ng = evaluate(d, run_dir)
    proc = subprocess.run(["Rscript", str(work / "edger.R"), str(run_dir)], capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(proc.stderr[-1500:])
    er = json.loads((run_dir / "edger.json").read_text())
    for kind in ("all", "support"):
        pv = np.full(len(truth), np.nan)
        for name, p in zip(er[kind]["region"], er[kind]["p"], strict=True):
            pv[int(name[1:])] = p
        res[f"C1_glm_{kind}"] = dict(p=pv.tolist(), estimate=[np.nan] * len(truth))
    print(f"  done {sc} {rep}", flush=True)
    return sc, rep, res, truth, ng


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--templates")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    templates = np.load(args.templates)["templates"] if args.templates else None
    settings = SimSettings(n_regions=300)
    work = Path(tempfile.mkdtemp(prefix="rm_", dir=out))
    (work / "edger.R").write_text(R_SCRIPT)
    jobs = [(sc, rep) for sc in SCENARIOS for rep in range(args.reps)]

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(run_one, [(sc, rep, work, settings, templates) for sc, rep in jobs]))
    summary = {}
    for sc, rep, res, truth, ng in results:
        kind = TRUTH[sc]
        for m, v in res.items():
            tr = np.asarray(truth, bool)
            if m.startswith("D_"):
                tr = np.array([tr[g * GROUP:(g + 1) * GROUP].any() for g in range(ng)])
            s = score(v["p"], tr)
            est = np.asarray(v["estimate"], float)
            s["median_estimate_changed"] = float(np.nanmedian(est[tr])) if tr.any() and np.isfinite(est).any() else float("nan")
            summary.setdefault(sc, {}).setdefault(m, []).append(s)
    table = {sc: {m: {k: float(np.nanmean([r[k] for r in rows])) for k in rows[0]} for m, rows in ms.items()}
             for sc, ms in summary.items()}
    (out / "region_methods.json").write_text(json.dumps(table, indent=1))
    lines = [f"# Per-region shape methods ({args.reps} runs, {'templates' if templates is not None else 'parametric'})", ""]
    methods = sorted(next(iter(table.values())).keys())
    lines += ["| method | " + " | ".join(SCENARIOS) + " |", "|---|" + "---|" * len(SCENARIOS)]
    for m in methods:
        cells = []
        for sc in SCENARIOS:
            t = table[sc][m]
            if TRUTH[sc] is None:
                cells.append(f"{t['called']:.3f} (any {t['any']:.2f})")
            else:
                cells.append(f"pow {t['power']:.2f} / FDP {t['fdp']:.2f} / est {t['median_estimate_changed']:.3g}")
        lines.append(f"| {m} | " + " | ".join(cells) + " |")
    (out / "REGION_METHODS.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
