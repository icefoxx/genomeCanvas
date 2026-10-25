from pathlib import Path

import numpy as np

from genomecanvas.io.bed import read_bed
from genomecanvas.io.cytoband import read_cytobands
from genomecanvas.io.gtf import read_gtf
from genomecanvas.layout.labels import label_layout
from genomecanvas.model.background_spec import BackgroundSpec
from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.model.scale import ScaleSpec
from genomecanvas.model.track import SignalSlice, Track
from genomecanvas.tracks.derived_signal import derive_signal
from genomecanvas.tracks.difference import compare_signal
from genomecanvas.transform.automatic_focus import largest_peak
from genomecanvas.transform.background import (
    estimate_background,
    rolling_background,
    transform_values,
)
from genomecanvas.validation.background_validator import validate_background
from genomecanvas.validation.figure_linter import lint_figure

FIXTURE = Path(__file__).parent / "fixtures/mini_genome"


def test_read_bed_skips_track_browser_headers(tmp_path: Path) -> None:
    path = tmp_path / "track.bed"
    path.write_text(
        "track type=bed name='HoxD'\n"
        "browser position chr2:1-10\n"
        "chr2\t10\t20\tpeak\t0\t+\n",
        encoding="utf-8",
    )
    features = read_bed(path)
    assert len(features) == 1
    assert features[0].interval.chrom == "chr2"


def test_fixture_readers() -> None:
    assert len(read_bed(FIXTURE / "peaks.bed")) == 2
    transcripts = read_gtf(FIXTURE / "mini.gtf", chrom="chr1")
    assert len(transcripts) == 2
    assert len(transcripts[0].exons) == 2


def test_difference_requires_explicit_ratio_pseudocount() -> None:
    starts = np.array([0, 100])
    ends = np.array([100, 200])
    left = SignalSlice("chr1", 0, 200, starts, ends, np.array([1.0, 2.0]), "counts")
    right = SignalSlice("chr1", 0, 200, starts, ends, np.array([2.0, 4.0]), "counts")
    assert np.allclose(compare_signal(left, right).values, [1.0, 2.0])


def test_compare_signal_accepts_chromosome_aliases() -> None:
    starts = np.array([0, 100])
    ends = np.array([100, 200])
    left = SignalSlice("chrX", 0, 200, starts, ends, np.array([1.0, 2.0]))
    right = SignalSlice("X", 0, 200, starts, ends, np.array([3.0, 5.0]))
    assert np.allclose(compare_signal(left, right).values, [2.0, 3.0])


def test_scale_modes_and_cytoband_coordinates() -> None:
    assert ScaleSpec("fixed", 0, 400).resolve(np.array([1.0, 2.0])) == (0, 400)
    bands = read_cytobands(
        Path(__file__).parent / "fixtures/reference_demo/mm10_cytoBandIdeo.txt.gz", chrom="chr7"
    )
    assert bands[0].start == 0
    assert bands[-1].end > bands[0].end


def test_local_background_estimators_and_transformations() -> None:
    left = np.full(10, 10.0)
    right = np.full(10, 10.0)
    right[0] = 110.0
    spec = BackgroundSpec(upstream=100, downstream=100)
    estimate = estimate_background(left, right, spec)
    assert estimate.value == 10.0
    assert estimate.estimator == "median"
    signal = np.array([110.0])
    assert np.allclose(transform_values(signal, estimate.value, "subtract"), [100.0])
    assert np.allclose(transform_values(signal, estimate.value, "ratio", 0.0), [11.0])
    assert np.allclose(transform_values(signal, estimate.value, "log2_ratio", 0.0), [np.log2(11)])


def test_derived_track_records_background_metadata() -> None:
    starts = np.array([100])
    ends = np.array([101])
    target = SignalSlice("chr1", 100, 101, starts, ends, np.array([110.0]), "counts")
    flank = SignalSlice("chr1", 0, 10, np.arange(10), np.arange(1, 11), np.full(10, 10.0), "counts")
    parent = Track("raw", "signal", "Raw signal", units="counts")
    derived = derive_signal(
        parent, target, flank, flank,
        BackgroundSpec(upstream=100, downstream=100), "log2_ratio", 0.01,
    )
    assert derived.metadata["background_estimator"] == "median"
    assert derived.metadata["transformation"] == "log2_ratio"


def test_background_validation_rejects_empty_flanks() -> None:
    issues = validate_background(BackgroundSpec(upstream=0, downstream=0))
    assert any(issue.severity == "ERROR" for issue in issues)


def test_advanced_scaling_focus_lint_and_labels() -> None:
    assert ScaleSpec("mad").resolve(np.array([1.0, 1.0, 1.0, 10.0]))[1] >= 1.0
    rolling = rolling_background(np.array([10.0, 110.0, 10.0]), 3)
    assert rolling[1] == 10.0
    starts = np.array([0, 100, 200])
    ends = np.array([100, 200, 300])
    signal = SignalSlice("chr1", 0, 300, starts, ends, np.array([1.0, 8.0, 2.0]), None)
    assert largest_peak(signal, 100).start == 100
    layout = label_layout([(0, 20, "high", 10), (10, 30, "low", 1)], leader_lines=True)
    assert len(layout) == 1 and layout[0]["leader_line"] is True


def test_figure_linter_warns_about_independent_comparison_scales() -> None:
    spec = FigureSpec.parse_obj({
        "genome": {"assembly": "hg38"},
        "region": {"chrom": "chr1", "start": 0, "end": 10},
        "tracks": [
            {"id": "a", "type": "signal", "source": {"path": "a.bg"}},
            {"id": "b", "type": "signal", "source": {"path": "b.bg"}},
        ],
        "scales": {"mode": "per_track"},
    })
    assert any(issue.severity == "WARNING" for issue in lint_figure(spec))
