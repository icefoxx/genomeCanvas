"""CLI for the deeptools-free heatmap pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..pipeline.run import HeatmapConfig, run_heatmap_pipeline


def run_heatmap(argv: list[str] | None = None) -> int:
    """Parse wrapper-compatible flags and run the pipeline."""
    parser = argparse.ArgumentParser(prog="gcanvas heatmap")
    parser.add_argument("-p", "--analysisPrefix", required=True, dest="prefix")
    parser.add_argument("-m", "--metaDataFile", required=True, dest="metadata")
    parser.add_argument("-r", "--regions", default="", dest="regions")
    parser.add_argument("-l", "--regionsLabels", default="auto", dest="region_labels")
    parser.add_argument("-f", "--flank", type=int, default=2000)
    parser.add_argument("--binsize", type=int, default=50)
    parser.add_argument("-s", "--start_point", choices=["center", "TSS", "TES"], default="center")
    parser.add_argument("--bodyLength", type=int, default=0, dest="body_length")
    parser.add_argument("--sortRegions", default="descend", dest="sort_regions")
    parser.add_argument("--sortUsing", default="mean", dest="sort_using")
    parser.add_argument("--sortUsingSamples", default="", dest="sort_using_samples")
    parser.add_argument("--extras", default=True, dest="extras")
    parser.add_argument("--overwrite", default=False, dest="overwrite")
    parser.add_argument(
        "--inputType", choices=["region", "bin"], default="region",
        dest="input_type",
    )
    parser.add_argument("--gBinSize", type=int, default=1000, dest="g_bin_size")
    parser.add_argument("--clusterRegions", default="none", dest="cluster_regions")
    parser.add_argument("--clusterNumber", type=int, default=5, dest="cluster_number")
    parser.add_argument("--clusterUsingSamples", default="", dest="cluster_using_samples")
    parser.add_argument("--clusteringDimReduct", default="raw", dest="clustering_dim")
    parser.add_argument(
        "--leidenResolutionParameter", type=float, default=1.0,
        dest="leiden_resolution",
    )
    parser.add_argument("--renormApproach", default="profileMinMax", dest="renorm_approach")
    parser.add_argument("--reuseRenormedMatrix", default=True, dest="reuse_renorm")
    parser.add_argument("--bigwig-engine", default="auto", dest="engine")
    parser.add_argument("--heatmapHeightBase", type=float, default=8.0, dest="heatmap_height")
    parser.add_argument("--rndSeed", default="7201", dest="seed")
    parser.add_argument("--diffPeakThresholds", default="Up2_Down2", dest="voom_threshold")
    parser.add_argument("--peakCenterBinFlank", type=int, default=5, dest="peak_center_flank")
    parser.add_argument(
        "--peakCenterSignifThreshold", type=float, default=0.05,
        dest="peak_center_alpha",
    )
    parser.add_argument("--maxRegionsPG", type=int, default=5000, dest="max_pairgrid_rows")
    parser.add_argument("--xAxisLabel", default="peak center", dest="x_label")
    parser.add_argument("--yAxisLabel", default="regions", dest="y_label")
    parser.add_argument("--title", default=None, dest="title")
    parser.add_argument("--colorbarLabel", default="signal", dest="colorbar_label")
    parser.add_argument("--sampleLabels", default="", dest="sample_labels")
    parser.add_argument(
        "--plotScale", choices=["linear", "log2", "zscore"], default="linear", dest="plot_scale",
    )
    parser.add_argument("--errorEstimates", default=True, dest="error_estimates")
    parser.add_argument("--lnV", default=True, dest="ln_v")
    parser.add_argument("--zMin", type=float, default=None, dest="zmin")
    parser.add_argument("--zMax", type=float, default=None, dest="zmax")
    args = parser.parse_args(argv)
    regions = _split_list(args.regions)
    labels = None if args.region_labels == "auto" else _split_list(args.region_labels)
    seed: int | None
    seed = None if str(args.seed) == "auto" else int(args.seed)
    config = HeatmapConfig(
        prefix=args.prefix,
        metadata=args.metadata,
        regions=regions,
        region_labels=labels,
        flank=args.flank,
        binsize=args.binsize,
        ref_point=args.start_point,
        body_length=args.body_length,
        sort_regions=args.sort_regions,
        sort_using=args.sort_using,
        sort_using_samples=args.sort_using_samples,
        extras=_as_bool(args.extras),
        overwrite=_as_bool(args.overwrite),
        input_type=args.input_type,
        g_bin_size=args.g_bin_size,
        cluster_regions=args.cluster_regions,
        cluster_number=args.cluster_number,
        cluster_using_samples=args.cluster_using_samples,
        clustering_dim=args.clustering_dim,
        leiden_resolution=args.leiden_resolution,
        renorm_approach=args.renorm_approach,
        reuse_renorm=_as_bool(args.reuse_renorm),
        engine=args.engine,
        heatmap_height=args.heatmap_height,
        seed=seed,
        voom_threshold=args.voom_threshold,
        peak_center_flank=args.peak_center_flank,
        peak_center_alpha=args.peak_center_alpha,
        max_pairgrid_rows=args.max_pairgrid_rows,
        x_label=args.x_label,
        y_label=args.y_label,
        title=args.title,
        colorbar_label=args.colorbar_label,
        sample_labels=_split_list(args.sample_labels) or None,
        plot_scale=args.plot_scale,
        error_estimates=_as_bool(args.error_estimates),
        ln_v=_as_bool(args.ln_v),
        zmin=args.zmin,
        zmax=args.zmax,
    )
    outputs = run_heatmap_pipeline(config)
    for name, path in outputs.items():
        print(f"{name}\t{path}")
    return 0


def _split_list(value: str) -> list[str]:
    if not value:
        return []
    if value.endswith(".lst") and Path(value).is_file():
        return [
            line.strip()
            for line in Path(value).read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    return [item.strip() for item in value.split(",") if item.strip()]


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() in {"1", "true", "yes", "t", "y"}
