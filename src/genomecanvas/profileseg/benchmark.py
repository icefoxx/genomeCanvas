"""Benchmark with known truth: false-positive rate, power, direction, bias and CI coverage.

The true effect of every metric in a scenario is obtained by running the identical analysis on
the noise-free expected profiles of that scenario, so bias and coverage are measured against
what the estimator targets, not against a hand-derived number.
"""

from __future__ import annotations

import time
from dataclasses import asdict

import numpy as np

from . import metrics as M
from .analysis import analyse
from .data import Config
from .simulate import SCENARIOS, SimSettings, simulate

METRICS = (
    "log2_auc",
    "log2_centre_auc",
    "log2_enrichment",
    "centroid",
    "sd_width",
    "iqr_width",
    "asymmetry",
    "naive_sd_width",
    "naive_centroid",
)
# A metric "truly changes" when its noise-free effect exceeds these tolerances. Consequences such
# as a translation lowering the ±250 bp centre AUC are real changes, not false positives.
TOLERANCE = dict(
    log2_enrichment=0.03,
    log2_auc=0.03,
    log2_centre_auc=0.03,
    centroid=3.0,
    sd_width=3.0,
    iqr_width=3.0,
    asymmetry=0.01,
    naive_sd_width=3.0,
    naive_centroid=3.0,
)


def true_effects(
    scenario: str, settings: SimSettings, config: Config, templates=None, seed: int = 0
) -> dict:
    profiles, _ = simulate(scenario, settings, seed=seed, templates=templates, noiseless=True)
    a = analyse(profiles, config, bootstrap=20)
    return {e["metric"]: e["estimate"] for e in a.effects}


def run_benchmark(
    scenarios=None,
    reps: int = 30,
    settings: SimSettings | None = None,
    config: Config | None = None,
    templates=None,
    seed: int = 1000,
    bootstrap: int = 200,
    progress=None,
) -> dict:
    settings = settings or SimSettings()
    config = config or Config(reference="WT", alternative="KO")
    scenarios = list(scenarios or SCENARIOS)
    records = []
    truths = {}
    began = time.time()
    for s_index, scenario in enumerate(scenarios):
        per_run_truth = []
        for rep in range(reps):
            run_seed = seed + 7919 * s_index + rep + 1
            t0 = time.time()
            profiles, _ = simulate(scenario, settings, seed=run_seed, templates=templates)
            a = analyse(profiles, config, bootstrap=bootstrap)
            elapsed = time.time() - t0
            # Truth for *this* run's regions: the estimand depends on the region set drawn.
            truth_run = true_effects(scenario, settings, config, templates, seed=run_seed)
            per_run_truth.append(truth_run)
            for e in a.effects:
                if e["metric"] not in METRICS:
                    continue
                truth = truth_run.get(e["metric"], np.nan)
                expected = int(np.sign(truth)) if abs(truth) > TOLERANCE[e["metric"]] else 0
                # Raw interval decision (comparable across metrics) and the role/Holm-adjusted call.
                welch_flag = bool(np.isfinite(e["welch_low"]) and (e["welch_low"] > 0 or e["welch_high"] < 0))
                adjusted_flag = e["call"] in ("higher", "lower")
                boot_flag = e["boot_low"] > 0 or e["boot_high"] < 0
                records.append(
                    dict(
                        scenario=scenario,
                        rep=rep,
                        seed=run_seed,
                        metric=e["metric"],
                        expected=expected,
                        exact_null=bool(abs(truth) <= 0.01 * TOLERANCE[e["metric"]]),
                        estimate=e["estimate"],
                        truth=truth,
                        welch_flag=welch_flag,
                        adjusted_flag=adjusted_flag,
                        role=e.get("role", ""),
                        boot_flag=boot_flag,
                        direction_ok=bool(welch_flag and np.sign(e["estimate"]) == expected)
                        if expected
                        else None,
                        covered=bool(e["welch_low"] <= truth <= e["welch_high"])
                        if np.isfinite(e["welch_low"])
                        else None,
                        perm_p=e["perm_p"],
                        perm_floor=e["perm_floor"],
                        seconds=elapsed,
                    )
                )
            if progress:
                progress(scenario, rep)
        truths[scenario] = {
            m: float(np.mean([t.get(m, np.nan) for t in per_run_truth])) for m in METRICS
        }
    summary = summarise(records)
    return dict(
        records=records,
        summary=summary,
        familywise=familywise(records),
        truths=truths,
        settings=asdict(settings),
        reps=reps,
        bootstrap=bootstrap,
        seconds=round(time.time() - began, 1),
        seed=seed,
    )


def summarise(records: list[dict]) -> list[dict]:
    out = []
    keys = sorted(
        {(r["scenario"], r["metric"]) for r in records},
        key=lambda k: (list(SCENARIOS).index(k[0]), METRICS.index(k[1])),
    )
    for scenario, metric in keys:
        rows = [r for r in records if r["scenario"] == scenario and r["metric"] == metric]
        n = len(rows)
        expected = int(round(np.mean([r["expected"] for r in rows]))) if rows else 0
        est = np.array([r["estimate"] for r in rows])
        truths_ = np.array([r["truth"] for r in rows])
        truth = float(np.mean(truths_))
        covered = [r["covered"] for r in rows if r["covered"] is not None]
        out.append(
            dict(
                scenario=scenario,
                metric=metric,
                expected=expected,
                exact_null=bool(all(r.get("exact_null", False) for r in rows)),
                n=n,
                welch_flag_rate=float(np.mean([r["welch_flag"] for r in rows])),
                boot_flag_rate=float(np.mean([r["boot_flag"] for r in rows])),
                direction_accuracy=float(np.mean(dirs))
                if expected
                and (dirs := [r["direction_ok"] for r in rows if r["direction_ok"] is not None])
                else np.nan,
                truth=truth,
                mean_estimate=float(est.mean()),
                bias=float(np.mean(est - truths_)),
                sd_estimate=float(est.std(ddof=1)) if n > 1 else np.nan,
                ci_coverage=float(np.mean(covered)) if covered else np.nan,
                n_covered=len(covered),
                perm_floor=float(rows[0]["perm_floor"]),
                mean_seconds=float(np.mean([r["seconds"] for r in rows])),
            )
        )
    return out


def familywise(records: list[dict], exact: bool = False) -> list[dict]:
    """Per scenario: share of runs with ≥1 false *adjusted* call among metrics whose truth is ≈ 0.

    Uses the reported calls (primary metrics unadjusted, secondary Holm-adjusted, comparators
    excluded), i.e. what a user would actually read off a report.
    """
    out = []
    for scenario in dict.fromkeys(r["scenario"] for r in records):
        runs: dict = {}
        for r in records:
            if r["scenario"] != scenario or r.get("role") == "comparator" or r["metric"].startswith("naive"):
                continue
            entry = runs.setdefault(r["rep"], dict(false=False, raw=False, eligible=0))
            if (r.get("exact_null", False) if exact else r["expected"] == 0):
                entry["eligible"] += 1
                entry["false"] |= bool(r.get("adjusted_flag"))
                entry["raw"] |= bool(r["welch_flag"])
        eligible = [v for v in runs.values() if v["eligible"]]
        if eligible:
            out.append(dict(
                scenario=scenario, runs=len(eligible),
                familywise_false_call_rate=float(np.mean([v["false"] for v in eligible])),
                familywise_rate_unadjusted=float(np.mean([v["raw"] for v in eligible])),
                null_metrics_per_run=float(np.mean([v["eligible"] for v in eligible])),
            ))
    return out


# ---------------------------------------------------------------- v2: per-region FDR / power


REGION_TOLERANCE = dict(log2_auc=0.1, log2_signal=0.1, log2_signal_relative=0.1, log2_enrichment=0.1, centroid=10.0,
                        mad_width=10.0, shift=10.0, shape_kernel=0.5, log2_scale=0.05)


def region_truth(scenario, settings, config, templates, seed) -> dict:
    """Per-region true effects from the noise-free profiles of the same run.

    Returns metric → (region indices, true difference). For the relative metric the reference level
    is the median change of regions the scenario leaves unchanged (its estimand: change vs the bulk).
    """
    from .regional import amplitude_matrices

    profiles, truth = simulate(scenario, settings, seed=seed, templates=templates, noiseless=True)
    a = analyse(profiles, config, bootstrap=20)
    n_ref = len(a.ref_idx)
    comp, amp = amplitude_matrices(a)
    out = {}
    diff = lambda m: m[n_ref:].mean(axis=0) - m[:n_ref].mean(axis=0)  # noqa: E731
    for name, m in amp.items():
        out[name] = (comp, diff(m))
    unchanged = truth["regions"]["amp_rel"][comp] == 0
    d = out["log2_signal"][1]
    out["log2_signal_relative"] = (comp, d - np.median(d[unchanged]))
    # Full-window moments, as in per_region_tests (regions with positive mass in every sample).
    pos = comp[np.all(a.core_mass[:, comp] > 0, axis=0)]
    e = a.core_edges
    p_full = M.normalise(a.core_signal[:, pos] * M.widths(e))
    out["mad_width"] = (pos, diff(M.mean_abs_deviation(p_full, e)))
    if profiles.oriented:
        out["centroid"] = (pos, diff(M.centroid(p_full, e)))
    # Shape tests: truth from the scenario (translation per region; broadening of every region).
    shift_sign = np.asarray(truth["regions"]["shift"])[comp]
    broad = scenario in ("broaden", "fragment_length")
    out["shift"] = (comp, shift_sign * float(settings.shift_bp))
    out["shape_kernel"] = (comp, ((shift_sign != 0) | broad).astype(float))
    return out


REGIONAL_SCENARIOS = ("null", "amplitude", "global_scale", "shift", "broaden", "sparse_shift",
                      "opposite_shift", "site_amplitude", "global_with_site")


def run_regional_benchmark(scenarios=REGIONAL_SCENARIOS, reps: int = 20, settings: SimSettings | None = None,
                           config: Config | None = None, templates=None, seed: int = 5000, fdr: float = 0.05) -> dict:
    from .regional import per_region_tests

    settings = settings or SimSettings()
    config = config or Config(reference="WT", alternative="KO", region_fdr=fdr)
    records = []
    for s_index, scenario in enumerate(scenarios):
        for rep in range(reps):
            run_seed = seed + 7919 * s_index + rep + 1
            profiles, _ = simulate(scenario, settings, seed=run_seed, templates=templates)
            a = analyse(profiles, config, bootstrap=20)
            res = per_region_tests(a, fdr=fdr)
            true = region_truth(scenario, settings, config, templates, run_seed)
            for metric, t in res.tests.items():
                idx, effect = true[metric]
                lookup = dict(zip(idx.tolist(), effect, strict=True))
                eff = np.array([lookup.get(int(i), np.nan) for i in t["region_index"]])
                tr = np.where(np.abs(np.nan_to_num(eff)) > REGION_TOLERANCE[metric], np.sign(eff), 0).astype(int)
                called = np.isin(t["call"], ("Higher", "Lower", "Changed"))
                sign = np.where(t["call"] == "Higher", 1, np.where(t["call"] == "Lower", -1, 0))
                if metric == "shape_kernel":  # omnibus: no direction
                    sign = np.where(called, tr, 0)
                false = called & (tr == 0)
                wrong_dir = called & (tr != 0) & (sign != tr)
                true_pos = called & (tr != 0) & (sign == tr)
                records.append(dict(
                    scenario=scenario, rep=rep, metric=metric, tested=int(len(tr)), discoveries=int(called.sum()),
                    false=int(false.sum()), wrong_direction=int(wrong_dir.sum()), true_changes=int((tr != 0).sum()),
                    fdp=float((false.sum() + wrong_dir.sum()) / max(called.sum(), 1)),
                    power=float(true_pos.sum() / (tr != 0).sum()) if (tr != 0).any() else np.nan,
                ))
    summary = []
    for scenario in scenarios:
        for metric in dict.fromkeys(r["metric"] for r in records):
            rows = [r for r in records if r["scenario"] == scenario and r["metric"] == metric]
            if not rows:
                continue
            powers = [r["power"] for r in rows if np.isfinite(r["power"])]
            summary.append(dict(
                scenario=scenario, metric=metric, runs=len(rows), tested=float(np.mean([r["tested"] for r in rows])),
                true_changes=float(np.mean([r["true_changes"] for r in rows])),
                discoveries=float(np.mean([r["discoveries"] for r in rows])),
                fdr=float(np.mean([r["fdp"] for r in rows])),
                runs_with_any_false=float(np.mean([(r["false"] + r["wrong_direction"]) > 0 for r in rows])),
                power=float(np.mean(powers)) if powers else np.nan,
            ))
    return dict(records=records, summary=summary, reps=reps, fdr_target=fdr, settings=asdict(settings))
