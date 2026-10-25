from pathlib import Path

import numpy as np
import pytest

from genomecanvas.io.bigwig import (
    binned_bigwig,
    bwtool_available,
    pybigwig_available,
    resolve_engine,
    set_bigwig_engine,
    using_bigwig_engine,
)
from genomecanvas.profile import extract_profile

DEMO = Path(__file__).parent / "fixtures/trackplot_demo/bw/H1_k4me3.bw"
PEAKS = Path(__file__).parent / "fixtures/trackplot_demo/narrowpeak/H1_Oct4.bed"


def test_resolve_engine_rejects_unknown_name() -> None:
    with pytest.raises(ValueError, match="engine must be"):
        resolve_engine("deeptools")  # type: ignore[arg-type]


def test_using_bigwig_engine_restores_previous() -> None:
    set_bigwig_engine(None)
    with using_bigwig_engine("pybigwig"):
        assert resolve_engine("auto", batch=False) == "pybigwig"
    assert resolve_engine("auto", batch=False) in {"pybigwig", "pybigtools", "bwtool"}


def test_explicit_engine_overrides_process_default() -> None:
    if not pybigwig_available():
        pytest.skip("pyBigWig is not installed")
    with using_bigwig_engine("bwtool" if bwtool_available() else "pybigwig"):
        assert resolve_engine("pybigwig") == "pybigwig"


@pytest.mark.skipif(not DEMO.is_file() or not pybigwig_available(), reason="demo BigWig missing")
def test_pybigwig_locus_bins_are_finite() -> None:
    values = binned_bigwig(
        str(DEMO), "chr6", 31_125_776, 31_144_789, n_bins=40, engine="pybigwig"
    )
    assert values.shape == (40,)
    assert np.isfinite(values).all()


@pytest.mark.skipif(
    not DEMO.is_file() or not bwtool_available() or not pybigwig_available(),
    reason="need demo BigWig, pyBigWig, and bwtool",
)
def test_bwtool_and_pybigwig_locus_bins_agree() -> None:
    kwargs = {
        "path": str(DEMO), "chrom": "chr6",
        "start": 31_125_776, "end": 31_144_789, "n_bins": 40,
    }
    python = binned_bigwig(**kwargs, engine="pybigwig")
    native = binned_bigwig(**kwargs, engine="bwtool")
    assert python.shape == native.shape
    mask = np.isfinite(python) & np.isfinite(native)
    correlation = np.corrcoef(python[mask], native[mask])[0, 1]
    assert correlation > 0.95


@pytest.mark.skipif(
    not DEMO.is_file() or not PEAKS.is_file() or not bwtool_available() or not pybigwig_available(),
    reason="need demo peaks and both engines",
)
def test_extract_profile_engine_switch_matches_shape() -> None:
    sources = {"k4me3": str(DEMO)}
    python = extract_profile(
        sources, bed=PEAKS, start_from="center", up=500, down=500, binsize=50, engine="pybigwig"
    )
    native = extract_profile(
        sources, bed=PEAKS, start_from="center", up=500, down=500, binsize=50, engine="bwtool"
    )
    assert python[0].values.shape == native[0].values.shape
    left = python[0].values.ravel()
    right = native[0].values.ravel()
    if np.std(left) > 0 and np.std(right) > 0:
        assert np.corrcoef(left, right)[0, 1] > 0.85


@pytest.mark.skipif(not bwtool_available(), reason="bwtool is required")
def test_bwtool_matrix_accepts_float_bed_scores(tmp_path: Path) -> None:
    from genomecanvas.io.bigwig import read_bwtool_matrix

    demo = Path(__file__).parent / "fixtures/seqplots_demo/GSM1208360_chrI_100Kb_q5_sample.bw"
    if not demo.is_file():
        pytest.skip("SeqPlots example BigWig missing")
    bed = tmp_path / "float_score.bed"
    bed.write_text("chrI\t20000\t21000\tpeak\t8.22581\t*\n", encoding="utf-8")
    matrix = read_bwtool_matrix(str(demo), bed, up=100, down=100, binsize=50)
    assert matrix.shape[0] == 1
    assert np.isfinite(matrix).any()


def test_bwtool_float_treats_na_as_zero() -> None:
    from genomecanvas.io.bigwig import _bwtool_float

    assert _bwtool_float("NA") == 0.0
    assert _bwtool_float("nan") == 0.0
    assert _bwtool_float("12.5") == 12.5
    assert _bwtool_float(None) == 0.0


@pytest.mark.skipif(
    not (Path(__file__).parent / "fixtures/deeptools_encode/bw/pgt_chrx_2.5e6_3.5e6.bw").is_file()
    or not bwtool_available(),
    reason="need fly BigWig and bwtool",
)
def test_bwtool_summary_remaps_chrX_to_X() -> None:
    from genomecanvas.analysis import extract_summary

    bigwig = Path(__file__).parent / "fixtures/deeptools_encode/bw/pgt_chrx_2.5e6_3.5e6.bw"
    bed = Path(__file__).parent / "fixtures/deeptools_encode/bed/chrX_windows_200.bed"
    summary = extract_summary({"fly": str(bigwig)}, bed, engine="bwtool")
    assert summary.values.shape[0] > 0
    assert np.isfinite(summary.values).all()
    assert float(np.nansum(summary.values)) > 0
