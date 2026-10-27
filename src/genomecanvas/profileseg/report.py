"""Run the v1 pipeline and write reports, TSV results, figures and a run log."""

from __future__ import annotations

import json
import platform
import shutil
import time
from importlib import metadata
from pathlib import Path

import numpy as np

from .analysis import Analysis, analyse
from .data import Config, Profiles, extract, read_deeptools_matrix, read_manifest, read_regions
from .figures import benchmark_figure, make_figures


def _versions():
    out = {"python": platform.python_version()}
    for name in ("genomecanvas", "numpy", "scipy", "matplotlib", "pyBigWig"):
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pass
    return out


def _tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    cols = list(rows[0])
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(cols) + "\n")
        for r in rows:
            handle.write("\t".join(_fmt(r[c]) for c in cols) + "\n")


def _fmt(v):
    if isinstance(v, float):
        return "NA" if not np.isfinite(v) else f"{v:.6g}"
    return str(v)


def load_inputs(config: Config) -> tuple[Profiles, list[str], list[dict]]:
    samples, flags = read_manifest(config.path(config.manifest))
    if config.matrix:
        profiles = (
            read_deeptools_matrix(config.path(config.matrix), samples)
            if not config.matrix.endswith(".npz")
            else Profiles.load(config.path(config.matrix))
        )
    else:
        if not config.regions:
            raise ValueError("set 'regions' (BED) or 'matrix' in the config")
        regions = read_regions(config.path(config.regions))
        profiles = extract(config, samples, regions)
    return profiles, flags, samples


def run(config: Config, output: str | Path | None = None) -> dict:
    began = time.time()
    out = Path(output) if output else config.path(config.output)
    for sub in ("reports", "results/metrics", "figures", "metadata"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    profiles, manifest_flags, samples = load_inputs(config)
    shutil.copy(config.path(config.manifest), out / "metadata" / "sample_manifest.tsv")
    profiles.save(out / "results" / "profiles.npz")
    if config.align_anchors:
        from .shapeclass import align_profiles
        profiles, anchor_shift, _ = align_profiles(profiles, config.background_width, config.core_flank,
                                                   config.align_max_shift, seed=config.seed,
                                                   estimator=config.background_estimator)
    a = analyse(profiles, config)
    if config.align_anchors:
        moved = anchor_shift != 0
        a.flags.append(f"anchors re-centred on condition-blind shape templates (max ±{config.align_max_shift} bp): "
                       f"{int(moved.sum())} regions moved, median |shift| {np.median(np.abs(anchor_shift[moved])) if moved.any() else 0:.0f} bp")
    from .figures import regional_figure
    from .regional import per_region_tests

    regional = per_region_tests(a, fdr=config.region_fdr)
    write_tables(a, out)
    write_regional(a, regional, out)
    figures = make_figures(a, out / "figures") + regional_figure(a, regional, out / "figures")
    write_reports(a, out, manifest_flags, samples)
    log = dict(
        config={k: v for k, v in vars(config).items()},
        versions=_versions(),
        seconds=round(time.time() - began, 2),
        n_samples=len(samples),
        n_regions=len(profiles.regions),
        included=int(a.included.sum()),
        figures=[str(p.relative_to(out)) for p in figures],
        flags=manifest_flags + a.flags,
    )
    (out / "run_log.json").write_text(json.dumps(log, indent=2, default=str))
    return dict(analysis=a, regional=regional, output=out, log=log)


def write_tables(a: Analysis, out: Path) -> None:
    m = out / "results" / "metrics"
    samples = [a.profiles.samples[i] for i in np.concatenate([a.ref_idx, a.alt_idx])]
    regions = a.profiles.regions
    rows = []
    for s, sample in enumerate(samples):
        for r, region in enumerate(regions):
            rows.append(
                dict(
                    sample_id=sample["sample_id"],
                    condition=sample["condition"],
                    region=region["name"],
                    chrom=region["chrom"],
                    start=region["start"],
                    end=region["end"],
                    included_in_shape=bool(a.included[r]),
                    background=float(a.background[s, r]),
                    **{k: float(v[s, r]) for k, v in a.amplitude.items()},
                )
            )
    _tsv(m / "amplitude_per_profile.tsv", rows)
    _tsv(
        m / "shape_per_replicate.tsv",
        [
            dict(
                sample_id=sample["sample_id"],
                condition=sample["condition"],
                aggregate_negative_fraction=float(a.aggregate_negative[s]),
                **{k: float(v[s]) for k, v in a.replicate_summaries.items()},
            )
            for s, sample in enumerate(samples)
        ],
    )
    _tsv(m / "effects.tsv", a.effects)
    inc = [regions[i] for i in np.flatnonzero(a.included)]
    _tsv(
        m / "region_effects.tsv",
        [
            dict(
                region=reg["name"],
                chrom=reg["chrom"],
                start=reg["start"],
                end=reg["end"],
                **{k: float(v[i]) for k, v in a.region_effects.items()},
            )
            for i, reg in enumerate(inc)
        ],
    )
    e = a.curves["edges"]
    _tsv(
        m / "group_distributions.tsv",
        [
            dict(
                bin_start=float(e[i]),
                bin_end=float(e[i + 1]),
                p_reference=float(a.curves["p_ref"][i]),
                p_alternative=float(a.curves["p_alt"][i]),
                difference=float(a.curves["diff"][i]),
                diff_low=float(a.curves["diff_low"][i]),
                diff_high=float(a.curves["diff_high"][i]),
            )
            for i in range(len(e) - 1)
        ],
    )
    (m / "group_comparison.json").write_text(
        json.dumps(dict(group=a.group, exclusion=a.exclusion), indent=2)
    )


def write_reports(a: Analysis, out: Path, manifest_flags: list[str], samples: list[dict]) -> None:
    c = a.config
    r = out / "reports"
    cond = a.profiles.conditions
    n_ref, n_alt = len(a.ref_idx), len(a.alt_idx)
    floor = a.effects[0]["perm_floor"]
    lines = [
        "# 00 · Discovery and inputs",
        "",
        f"* Input: {a.profiles.provenance.get('source')} · {len(samples)} samples · {len(a.profiles.regions)} regions",
        f"* Contrast: **{c.alternative} vs {c.reference}** ({n_alt} vs {n_ref} biological replicates)",
        f"* Conditions present: {', '.join(f'{k} ({int((cond == k).sum())})' for k in sorted(set(cond)))}",
        f"* Anchor: `{c.anchor}` · orientation: **{'oriented' if a.profiles.oriented else 'unoriented'}** "
        f"(signed metrics {'reported' if a.profiles.oriented else 'not applicable'})",
        f"* Grid: ±{c.flank} bp, {c.binsize} bp bins (bp coordinates relative to the anchor; minus-strand rows mirrored)",
        f"* Genome builds in manifest: {sorted({s.get('genome_build', '') for s in samples} - {''}) or 'not recorded'}",
        "",
        "## Flags",
        "",
    ]
    lines += [f"* {f}" for f in manifest_flags + a.flags] or ["* none"]
    (r / "00_discovery.md").write_text("\n".join(lines) + "\n")

    methods = f"""# 01 · Methods and assumptions

**Estimands are kept separate.**
- **Amplitude**: AUC (signal × bp) over the window and over ±{c.center_half} bp.
- **Shape**: background-corrected aggregate distribution.

Shape is measured as follows:
1. **Background:** per profile, the {c.background_estimator} of the outer {c.background_width} bp on each side is
   subtracted. (The mean is unbiased for skewed count noise. The median sits below it and leaves a positive floor,
   which makes taller peaks look narrower.)
2. **Core window:** ±{c.core_flank} bp.
3. **Inclusion (condition-independent):** a region is kept when its pooled (all-sample mean) background-corrected
   mass is positive and at or above the {c.min_mass_quantile:.0%} quantile (threshold
   {a.exclusion["threshold"]:.4g}); {a.exclusion["included"]} of {a.exclusion["total"]} regions were kept.
4. **Aggregation:** each region is divided by its pooled mass (a condition-independent factor), then averaged per
   replicate.
5. **Signed moments:** centroid, SD (including within-bin variance), IQR width and asymmetry come from the
   **signed** aggregate. Negative noise is not clipped, because clipping recreates the signal-to-noise confound.
   JSD uses the clipped profile and is descriptive only.
6. **Group distance:** W1/W2 between condition-mean aggregate distributions, computed exactly on bp edges.

**Inference (global target; unit = biological replicate):**
- each replicate is summarised across regions, then the conditions are contrasted;
- primary interval: Welch t on the replicate summaries ({c.confidence:.0%});
- secondary interval: hierarchical bootstrap ({c.bootstrap} draws), resampling replicates within condition and
  regions shared across conditions;
- exact permutation of replicate labels ({a.effects[0]["n_assignments"]} assignments).
  **The smallest attainable two-sided p is {floor:.3g}.**

A naive comparator (no background correction, whole window) is reported only to show the confound.
Per-region effects are descriptive (heterogeneity), not tests.

## Results

**Observed** (estimates, {c.alternative} − {c.reference}):

Primary metrics (pre-specified) are judged on their own Welch interval; secondary metrics are
Holm-adjusted together; naive metrics are comparators only. Bootstrap intervals are descriptive when a
condition has ≤ 3 replicates.

| metric | role | estimate | unit | Welch CI | Holm p | bootstrap CI | permutation p (floor) | replicates | call |
|---|---|---|---|---|---|---|---|---|---|
"""
    for e in a.effects:
        ci = (
            f"[{e['welch_low']:.3g}, {e['welch_high']:.3g}]" if np.isfinite(e["welch_low"]) else "—"
        )
        holm_p = f"{e['p_holm']:.3g}" if np.isfinite(e.get("p_holm", np.nan)) else "—"
        methods += (
            f"| {e['metric']} | {e.get('role', '')} | {e['estimate']:.3g} | {e['unit']} | {ci} | {holm_p} | "
            f"[{e['boot_low']:.3g}, {e['boot_high']:.3g}] | {e['perm_p']:.3g} ({e['perm_floor']:.3g}) | "
            f"{e.get('replicate_consistency', '')} | {e['call']} |\n"
        )
    calls = [e for e in a.effects if e["call"] in ("higher", "lower")]
    methods += "\n**Inferred:** "
    methods += (
        "; ".join(f"{e['metric']} {e['call']} ({e['estimate']:.3g} {e['unit']})" for e in calls)
        if calls
        else "no metric's Welch CI excludes zero"
    )
    methods += (
        ". These are replicate-level intervals. With this design permutation p cannot go below "
        f"{floor:.3g}, so no claim rests on permutation significance.\n\n"
        "**Hypothesised:** none generated automatically; biological interpretation requires the analyst.\n"
    )
    (r / "01_methods_and_assumptions.md").write_text(methods)

    neg = ", ".join(f"{v:.1%}" for v in a.aggregate_negative)
    qc = ["# 02 · QC", "", "## Region inclusion (condition-independent)", ""]
    qc += [f"* {k.replace('_', ' ')}: {v}" for k, v in a.exclusion.items()]
    qc += ["", "## Signal", ""]
    amp = a.amplitude
    for i, idx in enumerate(np.concatenate([a.ref_idx, a.alt_idx])):
        s = a.profiles.samples[idx]
        qc.append(
            f"* {s['sample_id']} ({s['condition']}): median AUC {np.nanmedian(amp['auc'][i]):.4g}, missing "
            f"{np.nanmean(amp['missing_fraction'][i]):.2%}, zero bins {np.nanmean(amp['zero_fraction'][i]):.2%}, "
            f"negative bins {np.nanmean(amp['negative_fraction'][i]):.2%}, "
            f"centre/flank {np.nanmedian(amp['centre_flank_ratio'][i]):.3g}"
        )
    qc += [
        "",
        f"Aggregate negative mass after background correction (per replicate): {neg}",
        "",
        "## Diagnostics",
        "",
    ]
    dg = a.diagnostics
    qc.append(f"* Background level, {c.alternative} vs {c.reference}: log2 ratio {dg.get('background_log2_ratio', float('nan')):+.3f} "
              "(|ratio| > 0.25 suggests normalisation / IP-efficiency differences; compare log2_enrichment)")
    if "amplitude_heterogeneity_ratio" in dg:
        qc.append(f"* Region-specific amplitude change / within-condition noise: {dg['amplitude_heterogeneity_ratio']:.2f} "
                  "(> 1.5: regions change unequally; aggregate shape can be re-weighted)")
    for key, h in dg.get("heterogeneity", {}).items():
        qc.append(f"* Region-level |Δ{key}|: median {h['median_abs_between']:.3g} between conditions vs "
                  f"{h['median_abs_within_split']:.3g} within-condition replicate split (ratio {h['ratio']:.2f})")
    if "concordance" in dg:
        qc += ["", "Replicate concordance (Spearman on region-centred log AUC; reflects region-specific variation only):", ""]
        names = dg["concordance"]["samples"]
        qc.append("| | " + " | ".join(names) + " |")
        qc.append("|---" * (len(names) + 1) + "|")
        for n_, row in zip(names, dg["concordance"]["rho_residual"], strict=True):
            qc.append(f"| {n_} | " + " | ".join(f"{v:.2f}" for v in row) + " |")
    qc += ["", "## Region-level heterogeneity (descriptive)", ""]
    for k, v in a.region_effects.items():
        v = v[np.isfinite(v)]
        if v.size:
            q = np.quantile(v, [0.1, 0.5, 0.9])
            qc.append(f"* {k}: 10/50/90 % = {q[0]:.3g} / {q[1]:.3g} / {q[2]:.3g} (n = {v.size})")
    (r / "02_QC.md").write_text("\n".join(qc) + "\n")
    write_limitations(a, out)


def write_limitations(a: Analysis, out: Path) -> None:
    c = a.config
    floor = a.effects[0]["perm_floor"]
    text = f"""# 04 · Limitations

* **Replicates:** {len(a.alt_idx)} vs {len(a.ref_idx)}. The exact permutation p-value cannot go below {floor:.3g}.
  Welch intervals on so few replicate summaries assume approximate normality. The v1 pilot benchmark showed mild
  under-coverage for some metrics (see 03_benchmark.md when available).
* **Global target only:** v1 asks whether a metric changes on average over the region set. It does not say which
  regions change (per-region tests are planned for v2) or where within the profile (segmentation and local
  inference are planned for v2). Per-region effects are descriptive.
* **Technical confounds are not removable by the model:**
  - fragment-length or read-extension differences change apparent width;
  - anchors chosen from one condition create apparent shifts;
  - bigWig smoothing alters shape.
  Check `bigwig_params` and `fragment_length` in the manifest (flags in 00_discovery.md).
* **Background:** estimated from the outer {c.background_width} bp. A neighbouring peak in the flank inflates it.
  The {c.background_estimator} estimator was used.
* **Normalisation:** signal is analysed as delivered. Shape metrics are scale-free; amplitude contrasts depend on
  upstream normalisation.
* **Orientation:** signed metrics (centroid, asymmetry) are only meaningful for oriented anchors; this run is
  {"oriented" if a.profiles.oriented else "unoriented"}.
* **Not implemented in v1:** count-based (BAM) models, segmentation, IWTomics comparison, equivalence testing.
"""
    (out / "reports" / "04_limitations.md").write_text(text)


def write_benchmark(result: dict, out: Path) -> list[Path]:
    out = Path(out)
    (out / "reports").mkdir(parents=True, exist_ok=True)
    (out / "results" / "metrics").mkdir(parents=True, exist_ok=True)
    _tsv(out / "results" / "metrics" / "benchmark_runs.tsv", result["records"])
    _tsv(out / "results" / "metrics" / "benchmark_summary.tsv", result["summary"])
    figs = benchmark_figure(result["summary"], out / "figures")
    s = result["settings"]
    lines = [
        "# 03 · Benchmark (simulation with known truth)",
        "",
        f"* {result['reps']} repetitions per scenario, {s['n_ref']} vs {s['n_alt']} replicates, {s['n_regions']} regions, "
        f"±{s['flank']} bp / {s['binsize']} bp, bootstrap {result['bootstrap']}, base seed {result['seed']}, "
        f"runtime {result['seconds']} s.",
        "* Truth for every metric is the same analysis applied to the noise-free expected profiles. A metric *truly "
        "changes* when that noise-free effect exceeds a tolerance (log2 0.03, bp 3, asymmetry 0.01), so consequential "
        "changes (e.g. a shift lowering the centre AUC) are not miscounted as false positives.",
        "* *naive* metrics (no background correction) measure raw shape including background; they are reported to "
        "show the signal-to-noise confound, not as an alternative estimator.",
        "* *call rate*: fraction of runs whose Welch CI excludes 0. For metrics that do not truly change it is the "
        "false-positive rate. *coverage*: fraction of Welch CIs containing the truth.",
        "",
        "| scenario | metric | truly changes | n | call rate (Welch) | call rate (bootstrap) | direction ok | truth | "
        "mean estimate | bias | CI coverage |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in result["summary"]:
        d = "—" if not r["expected"] else f"{r['direction_accuracy']:.2f}"
        lines.append(
            f"| {r['scenario']} | {r['metric']} | {'yes' if r['expected'] else 'no'} | {r['n']} | "
            f"{r['welch_flag_rate']:.2f} | {r['boot_flag_rate']:.2f} | {d} | {r['truth']:.3g} | "
            f"{r['mean_estimate']:.3g} | {r['bias']:.3g} | {r['ci_coverage']:.2f} |"
        )
    lines += ["", "## Family-wise false calls (what a report would show)", "",
              "| scenario | runs | metrics with truth ≈ 0 per run | ≥ 1 false call (primary + Holm) | ≥ 1 false call (unadjusted) |",
              "|---|---|---|---|---|"]
    for f in result.get("familywise", []):
        lines.append(f"| {f['scenario']} | {f['runs']} | {f['null_metrics_per_run']:.1f} | "
                     f"{f['familywise_false_call_rate']:.2f} | {f['familywise_rate_unadjusted']:.2f} |")
    lines += [
        "",
        "Pilot scale: rates have a binomial standard error of about "
        f"{100 * np.sqrt(0.05 * 0.95 / max(result['reps'], 1)):.0f} percentage points near 5 %. "
        "Run more repetitions before drawing firm conclusions.",
    ]
    (out / "reports" / "03_benchmark.md").write_text("\n".join(lines) + "\n")
    return figs


def write_regional(a: Analysis, regional, out: Path) -> None:
    """region_tests.tsv (one row per region × metric), summary, and reports/05_regional.md."""
    from .regional import REGION_METRICS, summary_rows

    regions = a.profiles.regions
    rows = []
    for name, res in regional.tests.items():
        for k, ridx in enumerate(res["region_index"]):
            reg = regions[int(ridx)]
            rows.append(dict(region=reg["name"], chrom=reg["chrom"], start=reg["start"], end=reg["end"],
                             strand=reg.get("strand", "."), metric=name, estimate=float(res["estimate"][k]),
                             moderated_t=float(res["t"][k]) if "t" in res else float("nan"),
                             df=float(res["df"][k]) if "df" in res else float("nan"), p=float(res["p"][k]),
                             fdr=float(res["fdr"][k]), call=str(res["call"][k]),
                             **({"shift_bp_xcorr": float(res["shift_bp_xcorr"][k]), "shift_bp_fit": float(res["shift_bp_fit"][k]),
                                 "translation_r2": float(res["translation_r2"][k]),
                                 "translation_ambiguity": float(res["translation_ambiguity"][k]),
                                 "translation_consistent": bool(res["translation_consistent"][k]),
                                 "translation_flag": str(res["translation_flag"][k]),
                                 "translation_fit": str(res["translation_fit"][k])}
                                if "translation_consistent" in res else
                                dict(shift_bp_xcorr=float("nan"), shift_bp_fit=float("nan"), translation_r2=float("nan"),
                                     translation_ambiguity=float("nan"), translation_consistent="", translation_flag="",
                                     translation_fit=""))))
    m = out / "results" / "metrics"
    _tsv(m / "region_tests.tsv", rows)
    summary = summary_rows(regional)
    _tsv(m / "region_tests_summary.tsv", summary)
    _tsv(m / "robust_global_shape.tsv", regional.robust_global)
    c = a.config
    agg = {e["metric"]: e for e in a.effects}
    lines = [
        "# 05 · Per-region tests (v2)", "",
        f"Each region is tested with replicates as observations ({len(a.alt_idx)} {c.alternative} vs "
        f"{len(a.ref_idx)} {c.reference}) using empirical-Bayes moderated t statistics (Smyth 2004; numerically "
        f"identical to limma::eBayes). Benjamini–Hochberg FDR < {c.region_fdr} per metric. Calls: Higher / Lower / "
        "Not_significant (`shape_kernel`: Changed / Not_significant). `Not_significant` is **not** evidence of no "
        "change (equivalence testing is not implemented).",
        "", "| metric | meaning | regions tested | higher | lower | changed | prior df | median estimate |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for s in summary:
        lines.append(f"| {s['metric']} | {s['description']} | {s['tested']} | {s['higher']} | {s['lower']} | "
                     f"{s['changed']} | {s['prior_df']:.3g} | {s['median_estimate']:.3g} |")
    lines += ["", f"Shape tests use the {regional.shape_regions} regions with positive pooled background-corrected "
              "signal. `shape_kernel` is the most powerful shape test (it uses the whole profile); `shift` gives the "
              "direction and a cross-correlation shift (bp) only where the translation description is consistent (`translation_consistent`; otherwise Changed, with the reason in `translation_flag`); consistency is not proof of a translation (one-sided gains can pass) and estimates near the detection limit are biased upwards; `shift_bp_fit` keeps the earlier template fit; `shape_any` combines the two shape tests per region before one BH adjustment (one list instead of two). `centroid` / `mad_width` (full-window moments) are kept for continuity and have low "
              "power on broad profiles.", ""]
    lines += [f"* {n}" for n in regional.notes]
    if regional.depth:
        dep = regional.depth
        lines += ["", f"## Depth ({dep['source']})", "",
                  "| sample | median fragments per shape-tested region |", "|---|---|"]
        lines += [f"| {k} | {v:.0f} |" for k, v in dep["median_fragments"].items()]
        lines += ["", f"Share of regions with fewer than {dep['threshold']} fragments in some sample: {dep['low_depth_share']:.0%}."]
    if regional.classes:
        lines += ["", "## Shape classes (condition-blind, descriptive)", "",
                  "Classes from a shift-only k-mean alignment of the pooled profiles (all samples, labels unused). "
                  "They do not enter any test. Per class: share of regions called, and among shift calls the "
                  "direction balance and the median |shift|.", "",
                  "| class | regions | shape_kernel called | shift called | higher / lower | higher share | median |shift| of calls (bp) |",
                  "|---|---|---|---|---|---|---|"]
        for r in regional.classes:
            lines.append(f"| {r['shape_class']} | {r['regions']} | {r.get('shape_kernel_called', float('nan')):.2f} | "
                         f"{r.get('shift_called', float('nan')):.2f} | {r.get('shift_higher', '—')} / {r.get('shift_lower', '—')} | "
                         f"{r.get('shift_higher_share', float('nan')):.2f} | {r.get('median_abs_shift_called_bp', float('nan')):.3g} |")
        _tsv(m / "shape_classes.tsv", regional.classes)
    lines += ["", "## Reading the amplitude tests", "",
              "* `log2_auc` (window AUC) and `log2_signal` (background-corrected core signal) are absolute and inherit "
              "upstream normalisation: a global change in IP efficiency makes most regions 'Lower' or 'Higher'.",
              "* `log2_signal_relative` removes a per-sample offset estimated on presumed-unchanged regions, so it finds "
              "regions that change **more than the bulk**. It assumes most regions do not change. It is built on "
              "background-corrected signal because window AUC contains background, which makes weak regions look "
              "relatively 'up' whenever IP efficiency drops.",
              "* `log2_enrichment` compares centre signal with the region's own background, which is robust to global rescaling.",
              f"* Spearman correlation of relative change with signal strength: {regional.strength_trend:+.2f} "
              "(|ρ| > 0.2 with many calls: suspect non-linear normalisation differences).",
              "", "## Robust global shape (median of per-region values, per replicate)", "",
              "| metric | estimate | Welch CI | regions | aggregate counterpart |", "|---|---|---|---|---|"]
    for rg in regional.robust_global:
        base = rg["metric"].replace("median_region_", "")
        counterpart = agg.get(base, {})
        ci = f"[{rg['welch_low']:.3g}, {rg['welch_high']:.3g}]" if np.isfinite(rg["welch_low"]) else "—"
        cp = f"{counterpart.get('estimate', float('nan')):.3g} ({counterpart.get('call', '')})" if counterpart else "—"
        lines.append(f"| {rg['metric']} | {rg['estimate']:.3g} | {ci} | {rg['n_regions']} | {cp} |")
        robust_change = np.isfinite(rg["welch_low"]) and (rg["welch_low"] > 0 or rg["welch_high"] < 0)
        related = ("sd_width", "iqr_width") if base == "sd_width" else (base,)
        called = [m for m in related if agg.get(m, {}).get("call") in ("higher", "lower")]
        if called and not robust_change:
            a.flags.append(f"aggregate {', '.join(called)} change is not supported by the median of per-region {base} "
                           f"(robust estimator, CI [{rg['welch_low']:.3g}, {rg['welch_high']:.3g}]): the aggregate change "
                           "may come from re-weighting of regions")
    lines += ["", "The median-of-regions estimator is unaffected by regions changing amplitude unequally; when it disagrees "
              "with the aggregate shape metric, prefer it for statements about a typical region."]
    (out / "reports" / "05_regional.md").write_text("\n".join(lines) + "\n")


def write_regional_benchmark(result: dict, out: Path) -> None:
    """reports/03b_benchmark_regional.md + region_benchmark_summary.tsv."""
    m = out / "results" / "metrics"
    m.mkdir(parents=True, exist_ok=True)
    _tsv(m / "region_benchmark_summary.tsv", result["summary"])
    lines = [
        "# 03b · Per-region tests: FDR and power (simulation)", "",
        f"{result['reps']} runs per scenario; target FDR {result['fdr_target']}. Truth per region comes from the "
        "noise-free profiles of the same run (|true change| > 0.1 log2 or > 10 bp counts as changed). "
        "FDR = mean false-discovery proportion (a discovery in the wrong direction counts as false); "
        "'≥1 false' = share of runs with any false discovery; power = share of truly changed regions called in "
        "the right direction.", "",
        "| scenario | metric | regions tested | truly changed | discoveries | FDR | ≥1 false | power |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for s in result["summary"]:
        power = f"{s['power']:.2f}" if np.isfinite(s["power"]) else "—"
        lines.append(f"| {s['scenario']} | {s['metric']} | {s['tested']:.0f} | {s['true_changes']:.1f} | "
                     f"{s['discoveries']:.1f} | {s['fdr']:.3f} | {s['runs_with_any_false']:.2f} | {power} |")
    lines += ["", "With few runs, FDR estimates are noisy (one false discovery in a run with few discoveries gives FDP = 1)."]
    (out / "reports").mkdir(parents=True, exist_ok=True)
    (out / "reports" / "03b_benchmark_regional.md").write_text("\n".join(lines) + "\n")
