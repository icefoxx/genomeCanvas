"""`gcanvas-profileseg`: replicate-aware amplitude/shape comparison of aligned profiles (v1 global + v2 per-region)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def _templates(config_path: str, limit: int):
    """Strongest reference-condition profiles of a real dataset, lightly smoothed, as shape templates."""
    from . import metrics as M
    from .data import load_config
    from .report import load_inputs
    from .simulate import _smooth

    cfg = load_config(config_path)
    profiles, _, _ = load_inputs(cfg)
    ref = profiles.values[profiles.conditions == cfg.reference]
    mean = ref.mean(axis=0)
    mean = mean[np.all(np.isfinite(mean), axis=1)]
    # Rank by background-corrected core mass (raw totals favour high-background regions) and
    # smooth lightly so per-region noise is not treated as true shape.
    bg = M.background(mean, profiles.edges, cfg.background_width, cfg.background_estimator)
    core = M.core_mask(profiles.edges, cfg.core_flank)
    signal = np.clip(mean[:, core] - bg[:, None], 0, None).sum(axis=1)
    order = np.argsort(signal)[::-1][:limit]
    print(f"templates: {len(order)} real profiles from {cfg.reference} ({config_path})")
    return _smooth(mean[order], 2.0), cfg.flank, cfg.binsize


def _benchmark(args, out: Path, templates=None, flank=2000, binsize=20):
    from .benchmark import run_benchmark, run_regional_benchmark
    from .data import Config
    from .report import write_benchmark, write_regional_benchmark
    from .simulate import SimSettings

    settings = SimSettings(
        n_ref=args.n_ref, n_alt=args.n_alt, n_regions=args.n_regions, flank=flank, binsize=binsize
    )
    done = {"n": 0}

    def progress(scenario, rep):
        done["n"] += 1
        if rep == args.reps - 1:
            print(f"  benchmark: {scenario} done", flush=True)

    result = run_benchmark(
        reps=args.reps,
        settings=settings,
        templates=templates,
        bootstrap=args.bootstrap,
        config=Config(reference="WT", alternative="KO", flank=flank, binsize=binsize),
        progress=progress,
    )
    write_benchmark(result, out)
    if getattr(args, "regional", True):
        regional = run_regional_benchmark(
            reps=args.reps, settings=settings, templates=templates,
            config=Config(reference="WT", alternative="KO", flank=flank, binsize=binsize),
        )
        print("  benchmark: per-region tests done", flush=True)
        write_regional_benchmark(regional, out)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(prog="gcanvas-profileseg", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    init = sub.add_parser("init", help="write an example config.yaml")
    init.add_argument("--out", default="config.yaml")
    run = sub.add_parser("run", help="analyse data described by a config")
    run.add_argument("--config", required=True)
    run.add_argument("--out", help="output directory (default: config 'output')")
    ex = sub.add_parser(
        "example", help="synthetic bigWigs → full analysis (+ optional pilot benchmark)"
    )
    ex.add_argument("--out", default="ProfileSeg_example")
    ex.add_argument("--scenario", default="shift")
    ex.add_argument(
        "--benchmark-reps", type=int, default=0, help="also run the benchmark with N repetitions"
    )
    bench = sub.add_parser("benchmark", help="simulation benchmark with known truth")
    bench.add_argument("--out", default="ProfileSeg_benchmark")
    bench.add_argument("--reps", type=int, default=30)
    bench.add_argument(
        "--templates-config", help="use real profile shapes from this config's data as templates"
    )
    bench.add_argument("--templates-max", type=int, default=2000)
    bench.add_argument("--no-regional", dest="regional", action="store_false",
                       help="skip the per-region (v2) FDR/power benchmark")
    pw = sub.add_parser("power", help="power guidance: simulate replicates × regions × effect sizes")
    pw.add_argument("--out", default="ProfileSeg_power")
    pw.add_argument("--reps", type=int, default=20, help="simulated experiments per cell")
    pw.add_argument("--replicates", type=int, nargs="+", default=[2, 3, 4, 6])
    pw.add_argument("--regions", type=int, nargs="+", default=[100, 300, 1000])
    pw.add_argument("--shift-bp", type=float, nargs="*", default=[25, 50, 100])
    pw.add_argument("--amplitude-fold", type=float, nargs="*", default=[1.2, 1.5])
    pw.add_argument("--broaden-factor", type=float, nargs="*", default=[1.2, 1.5])
    pw.add_argument("--site-fold", type=float, nargs="*", default=[1.5, 2.5])
    pw.add_argument("--templates-config", help="use real profile shapes from this config's data")
    pw.add_argument("--templates-max", type=int, default=2000)
    for p in (bench, ex):
        p.add_argument("--n-ref", type=int, default=3)
        p.add_argument("--n-alt", type=int, default=3)
        p.add_argument("--n-regions", type=int, default=300)
        p.add_argument("--bootstrap", type=int, default=200)
    args = parser.parse_args(argv)

    if args.command == "init":
        from .data import EXAMPLE_CONFIG

        Path(args.out).write_text(EXAMPLE_CONFIG)
        print(f"wrote {args.out}")
        return 0
    if args.command == "run":
        from .data import load_config
        from .report import run as run_pipeline

        result = run_pipeline(load_config(args.config), args.out)
        _summary(result)
        return 0
    if args.command == "example":
        from .data import load_config
        from .report import run as run_pipeline
        from .simulate import SimSettings, write_example

        out = Path(args.out)
        truth = write_example(
            out,
            SimSettings(n_ref=args.n_ref, n_alt=args.n_alt, n_regions=args.n_regions),
            scenario=args.scenario,
        )
        print(
            f"synthetic data ({args.scenario}): {out}/bigwig, regions.bed, metadata/sample_manifest.tsv"
        )
        result = run_pipeline(load_config(out / "config.yaml"), out / "results")
        (out / "results" / "truth.json").write_text(json.dumps(truth, indent=2, default=lambda v: np.asarray(v).tolist()))
        _summary(result)
        if args.benchmark_reps:
            args.reps = args.benchmark_reps
            _benchmark(args, out / "results")
            print(f"benchmark report: {out / 'results' / 'reports' / '03_benchmark.md'}")
        return 0
    if args.command == "power":
        from .power import power_grid, write_power
        from .simulate import SimSettings

        templates, flank, binsize, note = None, 2000, 20, "Parametric Gaussian peaks."
        if args.templates_config:
            templates, flank, binsize = _templates(args.templates_config, args.templates_max)
            note = f"Shapes from the strongest reference profiles of {args.templates_config}."
        effects = {k: v for k, v in (("shift_bp", args.shift_bp), ("amplitude_fold", args.amplitude_fold),
                                     ("broaden_factor", args.broaden_factor), ("site_fold", args.site_fold)) if v}

        def progress(row):
            print(f"  {row['effect']}={row['size']:g} n={row['replicates']} regions={row['regions']}: "
                  f"aggregate {row['aggregate_rate']:.2f} per-region {row['region_rate']:.2f}", flush=True)

        rows = power_grid(effects, args.replicates, args.regions, reps=args.reps,
                          base=SimSettings(flank=flank, binsize=binsize), templates=templates, progress=progress)
        write_power(rows, Path(args.out), note)
        print(f"power report: {Path(args.out) / 'reports' / '06_power.md'}")
        return 0
    if args.command == "benchmark":
        templates, flank, binsize = None, 2000, 20
        if args.templates_config:
            templates, flank, binsize = _templates(args.templates_config, args.templates_max)
        _benchmark(args, Path(args.out), templates, flank, binsize)
        print(f"benchmark report: {Path(args.out) / 'reports' / '03_benchmark.md'}")
        return 0
    return 1


def _summary(result):
    a, out = result["analysis"], result["output"]
    print(f"regions included in shape analysis: {a.exclusion['included']}/{a.exclusion['total']}")
    for e in a.effects:
        ci = f"[{e['welch_low']:.3g}, {e['welch_high']:.3g}]" if np.isfinite(e["welch_low"]) else ""
        print(f"  {e['metric']:16s} {e['estimate']:+9.3g} {e['unit']:10s} {ci:24s} {e['call']}")
    print(f"permutation p floor for this design: {a.effects[0]['perm_floor']:.3g}")
    print(f"reports: {out / 'reports'} · figures: {out / 'figures'}")


if __name__ == "__main__":
    sys.exit(main())
