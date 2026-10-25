from pathlib import Path

import numpy as np

from genomecanvas.pipeline.cluster import cluster_labels
from genomecanvas.pipeline.matrix import compute_region_matrix, sort_matrix
from genomecanvas.pipeline.metadata import read_sample_metadata
from genomecanvas.pipeline.plotting import plot_enrichment_heatmap, plot_profile_lines
from genomecanvas.pipeline.renorm import renormalize
from genomecanvas.pipeline.run import HeatmapConfig, run_heatmap_pipeline
from genomecanvas.pipeline.voom import extract_voom_beds

FIXTURE = Path(__file__).parent / "fixtures/mini_genome"


def _write_meta(path: Path) -> Path:
    meta = path / "meta.tsv"
    a = FIXTURE / "signal_a.bedgraph"
    b = FIXTURE / "signal_b.bedgraph"
    meta.write_text(
        "bigWig\tbigWigLabel\tconditionGrp\tmetaGrp\tcolorMap\tprofileColor\n"
        f"{a}\tA\tctrl\tmark\tBlues\t#d35400\n"
        f"{b}\tB\tctrl\tmark\tOranges\t#2980b9\n",
        encoding="utf-8",
    )
    return meta


def test_metadata_reader(tmp_path: Path) -> None:
    meta = _write_meta(tmp_path)
    samples = read_sample_metadata(meta)
    assert len(samples) == 2
    assert samples[0].label == "A"
    assert samples[1].meta_group == "mark"


def test_region_matrix_renorm_and_plots(tmp_path: Path) -> None:
    meta = _write_meta(tmp_path)
    samples = read_sample_metadata(meta)
    matrix = compute_region_matrix(
        samples,
        [(str(FIXTURE / "peaks.bed"), "peaks")],
        flank=50,
        binsize=25,
        engine="auto",
        skip_zeros=False,
    )
    assert matrix.values.shape[0] == 2
    assert len(matrix.sample_labels) == 2
    sorted_matrix = sort_matrix(matrix, method="mean")
    assert sorted_matrix.values.shape == matrix.values.shape
    scaled = renormalize(sorted_matrix, approach="profileMinMax")
    peaks = [float(np.nanmax(curve)) for curve in scaled.mean_profiles().values()]
    assert max(peaks) <= 1.0001
    plot_enrichment_heatmap(scaled, tmp_path / "hm.png")
    plot_profile_lines(scaled, tmp_path / "pf.png")
    assert (tmp_path / "hm.png").is_file()
    assert (tmp_path / "pf.png").is_file()


def test_kmeans_clustering(tmp_path: Path) -> None:
    meta = _write_meta(tmp_path)
    samples = read_sample_metadata(meta)
    matrix = compute_region_matrix(
        samples, [(str(FIXTURE / "peaks.bed"), "peaks")], flank=50, binsize=25, skip_zeros=False
    )
    labels = cluster_labels(matrix, method="kmeansExt", n_clusters=2, seed=1)
    assert labels.size == matrix.values.shape[0]
    assert set(labels).issubset({1, 2})


def test_voom_bed_extraction(tmp_path: Path) -> None:
    voom = tmp_path / "diff.vout.anno.Ranks.tsv"
    lines = ["Region\tRegulation\tp.value\tq.value\tlog2FC"]
    for i in range(12):
        lines.append(f"chr1:{i*10}-{i*10+5}\tUp\t0.001\t0.01\t2.0")
    for i in range(12):
        lines.append(f"chr1:{200+i*10}-{205+i*10}\tDown\t0.001\t0.01\t-2.0")
    voom.write_text("\n".join(lines) + "\n", encoding="utf-8")
    beds = extract_voom_beds(voom, tmp_path / "voom", threshold="Up2_Down2")
    assert len(beds) == 2


def test_pipeline_end_to_end(tmp_path: Path) -> None:
    meta = _write_meta(tmp_path)
    outputs = run_heatmap_pipeline(
        HeatmapConfig(
            prefix=str(tmp_path / "run"),
            metadata=str(meta),
            regions=[str(FIXTURE / "peaks.bed")],
            flank=50,
            binsize=25,
            extras=True,
            overwrite=True,
            skip_zeros=False,
            cluster_regions="none",
        )
    )
    assert outputs["heatmap"].is_file()
    assert outputs["profile"].is_file()
    assert outputs["renorm_heatmap"].is_file()
