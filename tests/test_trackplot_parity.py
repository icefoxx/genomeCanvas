from pathlib import Path

import numpy as np

from genomecanvas.analysis import diffpeak, extract_summary, pca_scores
from genomecanvas.io.chromhmm import read_chromhmm
from genomecanvas.io.genes import region_from_gtf
from genomecanvas.io.signal import query_signal
from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.profile import extract_profile, profile_heatmap, profile_plot, profile_summarize
from genomecanvas.render.renderer import render
from genomecanvas.transform.summarize import summarize_signals

FIXTURE = Path(__file__).parent / "fixtures/mini_genome"


def test_gene_symbol_resolves_local_gtf() -> None:
    interval = region_from_gtf(FIXTURE / "mini.gtf", "g1", padding=10)
    assert interval.chrom == "chr1"
    assert interval.start == 9
    assert interval.end == 270


def test_query_signal_bins_bedgraph() -> None:
    slice_ = query_signal(FIXTURE / "signal_a.bedgraph", "chr1", 0, 300, pixels=3)
    assert slice_.values.shape == (3,)
    assert slice_.values[0] == 1.0


def test_chromhmm_reader() -> None:
    states = read_chromhmm(FIXTURE / "chromhmm.bed")
    assert len(states) == 3
    assert states[0].state == "1"


def test_render_trackplot_track_types(tmp_path: Path) -> None:
    spec = FigureSpec.parse_obj({
        "genome": {"assembly": "synthetic"},
        "region": {"gene": "g1", "padding": 0},
        "gene_model": str(FIXTURE / "mini.gtf"),
        "collapse_transcripts": True,
        "highlights": [{"start": 50, "end": 80, "color": "#ffc41a", "alpha": 0.3}],
        "tracks": [
            {
                "id": "a", "type": "signal", "label": "A",
                "source": {"path": str(FIXTURE / "signal_a.bedgraph")},
                "style": {"color": "#d35400"}, "visibility": False,
            },
            {
                "id": "b", "type": "signal", "label": "B",
                "source": {"path": str(FIXTURE / "signal_b.bedgraph")},
                "style": {"color": "#2980b9"}, "visibility": False,
            },
            {"id": "overlay", "type": "overlay", "members": ["a", "b"], "height": 1.4},
            {
                "id": "peaks", "type": "peak",
                "source": {"path": str(FIXTURE / "peaks.bed")},
            },
            {
                "id": "hmm", "type": "chromhmm",
                "source": {"path": str(FIXTURE / "chromhmm.bed")},
                "height": 0.4,
            },
            {"id": "genes", "type": "gene", "height": 0.8},
            {
                "id": "cyto", "type": "cytoband",
                "source": {"path": str(FIXTURE / "cytoband.txt")},
                "height": 0.35,
            },
        ],
    })
    output = tmp_path / "parity.pdf"
    render(spec, output)
    assert output.is_file()
    assert spec.resolved().region.chrom == "chr1"


def test_profile_extract_summarize_and_heatmap(tmp_path: Path) -> None:
    matrices = extract_profile(
        {"A": str(FIXTURE / "signal_a.bedgraph"), "B": str(FIXTURE / "signal_b.bedgraph")},
        bed=FIXTURE / "regions.bed",
        start_from="center",
        up=40,
        down=40,
        binsize=10,
        input_type="bedgraph",
    )
    assert len(matrices) == 2
    assert matrices[0].values.shape[0] == 2
    summary = profile_summarize(matrices, conditions={"A": "grp", "B": "grp"})
    assert "grp" in summary
    assert profile_plot(summary, up=40, down=40, output=tmp_path / "profile.png")
    assert profile_heatmap(matrices, output=tmp_path / "heatmap.png")


def test_summary_pca_and_diffpeak() -> None:
    summary = extract_summary(
        {"A": str(FIXTURE / "signal_a.bedgraph"), "B": str(FIXTURE / "signal_b.bedgraph")},
        FIXTURE / "peaks.bed",
        input_type="bedgraph",
    )
    scores = pca_scores(summary)
    assert scores.coordinates.shape == (2, 2)
    result = diffpeak(
        summary, {"A": "ctrl", "B": "treat"}, numerator="treat", denominator="ctrl"
    )
    assert result.p_value.size == 2
    assert result.fdr.size == 2


def test_condition_summarize_signals() -> None:
    left = query_signal(FIXTURE / "signal_a.bedgraph", "chr1", 0, 300, pixels=3)
    right = query_signal(FIXTURE / "signal_b.bedgraph", "chr1", 0, 300, pixels=3)
    collapsed = summarize_signals({"A": left, "B": right}, {"A": "g", "B": "g"}, stat="mean")
    assert np.allclose(collapsed["g"].values, (left.values + right.values) / 2)
