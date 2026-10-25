"""Orchestrate the deeptools-free enrichment heatmap pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from ..io.bigwig import Engine
from .bins import bins_from_first_sample, summarize_bins
from .cluster import cluster_labels, recompute_cluster_matrix, write_cluster_beds
from .compare import peak_center_comparisons
from .correlation import correlation_heatmap, pair_grid, region_means_table, sample_pca
from .matrix import (
    compute_region_matrix,
    load_matrix,
    save_matrix,
    sort_matrix,
    write_mean_profiles,
)
from .metadata import SampleMeta, read_sample_metadata
from .plot_style import HeatmapPlotStyle, ProfilePlotStyle
from .plotting import plot_enrichment_heatmap, plot_profile_lines
from .renorm import renormalize
from .voom import extract_voom_beds

RefPoint = Literal["center", "TSS", "TES"]


@dataclass
class HeatmapConfig:
    """User-facing pipeline options (wrapper-compatible names)."""

    prefix: str
    metadata: str
    regions: list[str]
    region_labels: list[str] | None = None
    flank: int = 2000
    binsize: int = 50
    n_proc: int = 1
    ref_point: RefPoint = "center"
    body_length: int = 0
    sort_regions: str = "descend"
    sort_using: str = "mean"
    sort_using_samples: str = ""
    extras: bool = True
    overwrite: bool = False
    input_type: str = "region"
    g_bin_size: int = 1000
    cluster_regions: str = "none"
    cluster_number: int = 5
    cluster_using_samples: str = ""
    clustering_dim: str = "raw"
    leiden_resolution: float = 1.0
    renorm_approach: str = "profileMinMax"
    reuse_renorm: bool = True
    engine: Engine = "auto"
    skip_zeros: bool = True
    zmin: float | None = None
    zmax: float | None = None
    heatmap_height: float = 8.0
    seed: int | None = 7201
    voom_threshold: str = "Up2_Down2"
    peak_center_flank: int = 5
    peak_center_alpha: float = 0.05
    max_pairgrid_rows: int = 5000
    x_label: str = "peak center"
    y_label: str = "regions"
    title: str | None = None
    colorbar_label: str = "signal"
    sample_labels: list[str] | None = None
    plot_scale: Literal["linear", "log2", "zscore"] = "linear"
    error_estimates: bool = True
    ln_v: bool = True


def run_heatmap_pipeline(config: HeatmapConfig) -> dict[str, Path]:
    """Run region or bin-mode enrichment plots. Returns written artifact paths."""
    samples = read_sample_metadata(config.metadata)
    prefix = Path(config.prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    if config.seed is not None:
        np.random.seed(config.seed)

    beds = _resolve_regions(config)
    if config.input_type == "bin":
        return _run_bin_mode(config, samples, prefix, outputs)

    matrix_path = Path(f"{prefix}.matrix.npz")
    if matrix_path.is_file() and not config.overwrite:
        matrix = load_matrix(matrix_path)
    else:
        matrix = compute_region_matrix(
            samples,
            beds,
            flank=config.flank,
            binsize=config.binsize,
            ref_point=config.ref_point,
            body_length=config.body_length if config.ref_point in {"TSS", "TES"} else 0,
            skip_zeros=config.skip_zeros,
            engine=config.engine,
        )
        if config.sort_regions in {"descend", "ascend"}:
            indices = _sample_indices(config.sort_using_samples, len(samples))
            matrix = sort_matrix(
                matrix,
                method=config.sort_using,  # type: ignore[arg-type]
                descending=config.sort_regions == "descend",
                sample_indices=indices,
            )
        save_matrix(matrix, matrix_path)
    outputs["matrix"] = matrix_path
    outputs["mean_profile"] = write_mean_profiles(
        matrix, Path(f"{prefix}.MeanSampleEnrichment.tsv")
    )
    heatmap_style, profile_style = _plot_styles(config)
    outputs["heatmap"] = plot_enrichment_heatmap(
        matrix, Path(f"{prefix}.heatmap.png"), style=heatmap_style,
    )
    outputs["profile"] = plot_profile_lines(
        matrix, Path(f"{prefix}.profile.png"), style=profile_style,
    )

    if config.extras and matrix.values.shape[0] > 10:
        table = region_means_table(matrix.values, matrix.sample_boundaries)
        outputs["pearson"] = correlation_heatmap(
            table, matrix.sample_labels, Path(f"{prefix}.corr.pearson.png"), method="pearson"
        )
        outputs["spearman"] = correlation_heatmap(
            table, matrix.sample_labels, Path(f"{prefix}.corr.spearman.png"), method="spearman"
        )
        outputs["pca"] = sample_pca(table, matrix.sample_labels, Path(f"{prefix}.pca.png"))
        outputs["pairgrid"] = pair_grid(
            table, matrix.sample_labels, Path(f"{prefix}.PairGrid.png"),
            max_rows=config.max_pairgrid_rows,
        )

    if len(matrix.group_labels) >= 2:
        written = peak_center_comparisons(
            matrix, prefix,
            flank_bins=config.peak_center_flank,
            significance=config.peak_center_alpha,
        )
        for path in written:
            outputs[path.name] = path

    renorm_path = Path(f"{prefix}.renorm.matrix.npz")
    if renorm_path.is_file() and config.reuse_renorm and not config.overwrite:
        renormed = load_matrix(renorm_path)
    else:
        renormed = renormalize(matrix, approach=config.renorm_approach)
        save_matrix(renormed, renorm_path)
    outputs["renorm_matrix"] = renorm_path
    outputs["renorm_heatmap"] = plot_enrichment_heatmap(
        renormed, Path(f"{prefix}.heatmap.renorm.png"),
        style=_plot_styles(config, zmin=0.0, zmax=1.1)[0],
    )
    outputs["renorm_heatmap_auto"] = plot_enrichment_heatmap(
        renormed, Path(f"{prefix}.heatmap.renorm.auto.png"),
        style=_plot_styles(config, zmin=0.0, autoscale=True)[0],
    )
    outputs["renorm_profile"] = plot_profile_lines(
        renormed, Path(f"{prefix}.profile.renorm.png"), style=profile_style,
    )

    if config.cluster_regions not in {"none", ""}:
        if len(beds) != 1:
            raise ValueError("clustering requires exactly one region BED")
        labels = cluster_labels(
            renormed,
            method=config.cluster_regions,  # type: ignore[arg-type]
            n_clusters=config.cluster_number,
            sample_indices=_sample_indices(config.cluster_using_samples, len(samples)),
            dim=config.clustering_dim,
            seed=config.seed,
            leiden_resolution=config.leiden_resolution,
        )
        cluster_prefix = Path(f"{prefix}.{config.cluster_regions}")
        cluster_beds = write_cluster_beds(renormed, labels, cluster_prefix)
        cluster_matrix = recompute_cluster_matrix(
            samples,
            cluster_beds,
            flank=config.flank,
            binsize=config.binsize,
            ref_point=config.ref_point,
            body_length=config.body_length if config.ref_point in {"TSS", "TES"} else 0,
            engine=config.engine,
        )
        cluster_renorm = renormalize(cluster_matrix, approach=config.renorm_approach)
        outputs["cluster_heatmap"] = plot_enrichment_heatmap(
            cluster_renorm, Path(f"{cluster_prefix}.heatmap.png"),
            style=_plot_styles(config, zmin=0.0, zmax=1.1)[0],
        )
        outputs["cluster_profile"] = plot_profile_lines(
            cluster_renorm, Path(f"{cluster_prefix}.profile.png"), style=profile_style,
        )
        save_matrix(cluster_renorm, Path(f"{cluster_prefix}.matrix.npz"))
    return outputs


def _plot_styles(
    config: HeatmapConfig,
    *,
    zmin: float | None = None,
    zmax: float | None = None,
    autoscale: bool | None = None,
) -> tuple[HeatmapPlotStyle, ProfilePlotStyle]:
    """Build heatmap/profile chrome from pipeline config."""
    vmin = config.zmin if zmin is None else zmin
    vmax = config.zmax if zmax is None else zmax
    auto = config.zmin is None and config.zmax is None
    if autoscale is not None:
        auto = autoscale
    heatmap = HeatmapPlotStyle(
        title=config.title,
        labels=config.sample_labels,
        xlab=config.x_label,
        ylab=config.y_label,
        colorbar_label=config.colorbar_label,
        zmin=vmin,
        zmax=vmax,
        autoscale=auto,
        height=config.heatmap_height,
        plot_scale=config.plot_scale,
        error_estimates=config.error_estimates,
        ln_v=config.ln_v,
    )
    profile = ProfilePlotStyle(
        title=config.title,
        labels=config.sample_labels,
        xlab=config.x_label,
        error_estimates=config.error_estimates,
        ln_v=config.ln_v,
        plot_scale=config.plot_scale,
    )
    return heatmap, profile


def _run_bin_mode(
    config: HeatmapConfig,
    samples: list[SampleMeta],
    prefix: Path,
    outputs: dict[str, Path],
) -> dict[str, Path]:
    bed = bins_from_first_sample(samples, config.g_bin_size, Path(f"{prefix}.bins.bed"))
    labels, table = summarize_bins(samples, bed, engine=config.engine)
    if config.extras and table.shape[0] > 10 and table.shape[1] > 1:
        outputs["pearson"] = correlation_heatmap(
            table, labels, Path(f"{prefix}.corr.pearson.png"), method="pearson"
        )
        outputs["spearman"] = correlation_heatmap(
            table, labels, Path(f"{prefix}.corr.spearman.png"), method="spearman"
        )
        outputs["pca"] = sample_pca(table, labels, Path(f"{prefix}.pca.png"))
        outputs["pairgrid"] = pair_grid(
            table, labels, Path(f"{prefix}.PairGrid.png"), max_rows=config.max_pairgrid_rows
        )
    outputs["bins"] = bed
    return outputs


def _resolve_regions(config: HeatmapConfig) -> list[tuple[str, str]]:
    regions = list(config.regions)
    labels = list(config.region_labels or [])
    if len(regions) == 1 and regions[0].endswith(".vout.anno.Ranks.tsv"):
        return extract_voom_beds(
            regions[0], config.prefix, threshold=config.voom_threshold, seed=config.seed or 7201
        )
    if len(regions) == 1 and regions[0].endswith(".lst"):
        regions = [
            line.strip() for line in Path(regions[0]).read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    if labels == ["auto"] or not labels:
        labels = [Path(path).stem for path in regions]
    if len(labels) != len(regions):
        raise ValueError("region labels must match the number of BED files")
    missing = [path for path in regions if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError("region files are not accessible: " + ", ".join(missing))
    return list(zip(regions, labels, strict=True))


def _sample_indices(spec: str, n_samples: int) -> list[int] | None:
    text = spec.strip()
    if not text:
        return None
    indices = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        value = int(part)
        if value < 1 or value > n_samples:
            raise ValueError(f"sample index {value} is outside 1..{n_samples}")
        indices.append(value - 1)
    return indices or None
