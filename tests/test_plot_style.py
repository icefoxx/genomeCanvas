from pathlib import Path

import numpy as np

from genomecanvas.pipeline.matrix import compute_region_matrix
from genomecanvas.pipeline.metadata import SampleMeta
from genomecanvas.pipeline.plot_style import (
    HeatmapPlotStyle,
    ProfilePlotStyle,
    format_distance,
    format_genomic_position,
    scale_values,
)
from genomecanvas.pipeline.plotting import plot_enrichment_heatmap, plot_profile_lines

FIXTURE = Path(__file__).parent / "fixtures/mini_genome"


def test_format_distance_auto_units() -> None:
    assert format_distance(0) == "0"
    assert "kb" in format_distance(2000)
    assert format_distance(2_000_000, "Mb").endswith("Mb")


def test_format_genomic_position_uses_commas_for_bp() -> None:
    assert "," in format_genomic_position(31_125_776, span=2000, unit="bp")


def test_scale_values_log2_and_zscore() -> None:
    values = np.array([[0.0, 3.0], [1.0, 7.0]])
    logged = scale_values(values, "log2")
    assert np.isfinite(logged).all()
    zscored = scale_values(values, "zscore")
    assert abs(float(np.nanmean(zscored))) < 1e-6


def test_heatmap_and_profile_accept_custom_labels(tmp_path: Path) -> None:
    samples = [
        SampleMeta(str(FIXTURE / "signal_a.bedgraph"), "A", "g", "g", "Blues", "auto"),
        SampleMeta(str(FIXTURE / "signal_b.bedgraph"), "B", "g", "g", "Oranges", "#d35400"),
    ]
    matrix = compute_region_matrix(
        samples, [(str(FIXTURE / "peaks.bed"), "peaks")],
        flank=50, binsize=25, skip_zeros=False,
    )
    plot_enrichment_heatmap(
        matrix,
        tmp_path / "hm.png",
        style=HeatmapPlotStyle(
            title="Custom title",
            labels=("Control", "Treat"),
            xlab="from peak (bp)",
            ylab="peaks",
            colorbar_label="coverage",
            height=5.0,
        ),
    )
    plot_profile_lines(
        matrix,
        tmp_path / "pf.png",
        style=ProfilePlotStyle(
            title="Mean coverage",
            labels=("Control", "Treat"),
            xlab="from peak (bp)",
            ylab="mean",
        ),
    )
    assert (tmp_path / "hm.png").is_file()
    assert (tmp_path / "pf.png").is_file()
