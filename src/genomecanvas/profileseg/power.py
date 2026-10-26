"""Power guidance: simulate a grid of designs and report how often each effect is detected.

Each cell (replicates per condition × number of regions × effect) is simulated ``reps`` times with
the same noise model as the benchmark, optionally on real profile shapes (templates). For every
cell the matching false-call rate (same design, no effect) is reported alongside power, so a
design is only "adequate" when power is high *and* false calls stay near the nominal level.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from .analysis import analyse
from .data import Config
from .regional import per_region_tests
from .simulate import SimSettings, simulate

# effect name → (scenario, SimSettings field, aggregate metric, per-region metric or None)
EFFECTS = {
    "shift_bp": ("shift", "shift_bp", "centroid", "centroid"),
    "amplitude_fold": ("amplitude", "amplitude", "log2_auc", None),
    "broaden_factor": ("broaden", "broaden", "iqr_width", None),
    "site_fold": ("site_amplitude", "site_fold", None, "log2_signal_relative"),
}
NULL_SETTING = {"shift_bp": 0.0, "amplitude_fold": 1.0, "broaden_factor": 1.0, "site_fold": 1.0}


def _calls(a, metric):
    e = next((x for x in a.effects if x["metric"] == metric), None)
    if e is None or not np.isfinite(e["welch_low"]):
        return False, np.nan
    return bool(e["welch_low"] > 0 or e["welch_high"] < 0), float(e["estimate"])


def power_cell(effect: str, size: float, n: int, regions: int, reps: int, base: SimSettings,
               templates=None, seed: int = 0, fdr: float = 0.05) -> dict:
    scenario, field_name, agg_metric, reg_metric = EFFECTS[effect]
    settings = replace(base, n_ref=n, n_alt=n, n_regions=regions, **{field_name: size})
    config = Config(reference="WT", alternative="KO", flank=base.flank, binsize=base.binsize, region_fdr=fdr)
    agg_hits, reg_power, reg_fdp = [], [], []
    for rep in range(reps):
        profiles, truth = simulate(scenario, settings, seed=seed + rep, templates=templates)
        a = analyse(profiles, config, bootstrap=10)
        if agg_metric:
            agg_hits.append(_calls(a, agg_metric)[0])
        if reg_metric:
            t = per_region_tests(a, fdr=fdr).tests.get(reg_metric)
            if t is None:
                continue
            key = "shift" if reg_metric == "centroid" else "amp_rel"
            true = truth["regions"][key][t["region_index"]]
            if size == NULL_SETTING[effect]:
                true = np.zeros_like(true)  # the scenario's region mask is meaningless without an effect
            called = np.isin(t["call"], ("Higher", "Lower"))
            sign = np.where(t["call"] == "Higher", 1, -1)
            if (true != 0).any():
                reg_power.append(float((called & (sign == true) & (true != 0)).sum() / (true != 0).sum()))
            reg_fdp.append(float((called & (true == 0)).sum() / max(called.sum(), 1)))
    return dict(effect=effect, size=size, replicates=n, regions=regions, reps=reps,
                aggregate_metric=agg_metric or "", aggregate_rate=float(np.mean(agg_hits)) if agg_hits else np.nan,
                region_metric=reg_metric or "", region_rate=float(np.mean(reg_power)) if reg_power else np.nan,
                region_fdr=float(np.mean(reg_fdp)) if reg_fdp else np.nan)


def power_grid(effects: dict, replicates, regions, reps: int = 20, base: SimSettings | None = None,
               templates=None, seed: int = 31_000, progress=None) -> list[dict]:
    """``effects`` maps effect name → list of sizes. A null cell (no effect) is added for each design."""
    base = base or SimSettings()
    rows = []
    for e_index, (effect, sizes) in enumerate(effects.items()):
        for n in replicates:
            for m in regions:
                for size in [NULL_SETTING[effect], *sizes]:
                    rows.append(power_cell(effect, size, n, m, reps, base, templates,
                                           seed=seed + 100_000 * e_index + 1000 * n + m))
                    if progress:
                        progress(rows[-1])
    return rows


def write_power(rows: list[dict], out, templates_note: str = "") -> None:
    from pathlib import Path

    from .report import _tsv

    out = Path(out)
    (out / "results" / "metrics").mkdir(parents=True, exist_ok=True)
    (out / "reports").mkdir(parents=True, exist_ok=True)
    _tsv(out / "results" / "metrics" / "power_grid.tsv", rows)
    lines = ["# 06 · Power guidance (simulation)", "",
             f"{rows[0]['reps'] if rows else 0} simulated experiments per cell. {templates_note}",
             "Aggregate rate = share of experiments whose Welch interval excludes 0 for the named metric; at the null size "
             "this is the false-call rate. Per-region power = share of truly changed regions called (BH FDR 0.05); "
             "per-region FDR is shown alongside.", "",
             "| effect | size | replicates / condition | regions | aggregate metric | aggregate rate | per-region metric | per-region power | per-region FDR |",
             "|---|---|---|---|---|---|---|---|---|"]
    fmt = lambda v: "—" if not np.isfinite(v) else f"{v:.2f}"  # noqa: E731
    for r in rows:
        null = " (null)" if r["size"] == NULL_SETTING[r["effect"]] else ""
        lines.append(f"| {r['effect']} | {r['size']:g}{null} | {r['replicates']} | {r['regions']} | {r['aggregate_metric'] or '—'} | "
                     f"{fmt(r['aggregate_rate'])} | {r['region_metric'] or '—'} | {fmt(r['region_rate'])} | {fmt(r['region_fdr'])} |")
    lines += ["", "Rules of thumb are only valid for the noise level simulated here; run with `--templates-config` on a "
              "pilot of your own data to calibrate signal-to-noise."]
    (out / "reports" / "06_power.md").write_text("\n".join(lines) + "\n")
