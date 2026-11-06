"""Scale modes of enrichment heatmaps: colour (z) scale and profile y axis.

Both use scale_groups(): auto / shared / independent / group, plus per-sample user limits.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from genomecanvas.cli.heatmap import _float_list
from genomecanvas.pipeline import plotting
from genomecanvas.pipeline.matrix import EnrichmentMatrix, RegionCoord
from genomecanvas.pipeline.plot_style import HeatmapPlotStyle
from genomecanvas.pipeline.plotting import heatmap_limits, plot_enrichment_heatmap, scale_groups


def _matrix(peaks: list[float], n_rows: int = 20, n_bins: int = 10, groups=None) -> EnrichmentMatrix:
    shape = np.exp(-((np.arange(n_bins) - n_bins / 2) ** 2) / 4)
    blocks = [np.tile(1 + peak * shape, (n_rows, 1)) for peak in peaks]
    n = len(peaks)
    return EnrichmentMatrix(
        values=np.hstack(blocks),
        coords=[RegionCoord("chr1", i * 100, i * 100 + 10, ".", ".", ".", "g") for i in range(n_rows)],
        sample_labels=tuple(f"s{i}" for i in range(n)),
        sample_paths=tuple(f"s{i}.bw" for i in range(n)),
        sample_boundaries=tuple(range(0, n_bins * (n + 1), n_bins)),
        group_labels=("g",), group_boundaries=(0, n_rows),
        upstream=500, downstream=500, body=0, bin_size=100, ref_point="center",
        bins_per_sample=n_bins, meta_groups=tuple(groups or ["m"] * n), conditions=("c",) * n,
        colormaps=("Reds",) * n, profile_colors=("#aa0000",) * n,
    )


def _blocks(matrix: EnrichmentMatrix) -> list:
    return [matrix.sample_block(i) for i in range(len(matrix.sample_labels))]


def _draw(monkeypatch, matrix, style, tmp_path: Path):
    captured = {}
    original = plotting._apply_profile_ylims

    def keep(axes, opts, meta_groups=None):
        original(axes, opts, meta_groups)
        captured["limits"] = [ax.get_ylim() for ax in axes]
        captured["labels"] = [ax.yaxis.get_major_ticks()[0].label1.get_visible() for ax in axes]

    monkeypatch.setattr(plotting, "_apply_profile_ylims", keep)
    plot_enrichment_heatmap(matrix, tmp_path / "hm.png", style=style)
    plt.close("all")
    return captured


# --- grouping -------------------------------------------------------------------------------

def test_scale_groups_modes() -> None:
    values = [16.0, 18.0, 4.0, 8.0, 3.5]
    assert scale_groups(values, "shared") == [0, 0, 0, 0, 0]
    assert scale_groups(values, "independent") == [0, 1, 2, 3, 4]
    # auto, 2-fold: {3.5, 4} | {8} | {16, 18}; ids numbered in column order
    assert scale_groups(values, "auto", 2.0) == [0, 0, 1, 2, 1]
    # auto, 3-fold: cut at the largest fold gap (6 -> 16): {3.5, 4, 6} | {16, 18}
    assert scale_groups([16.0, 18.0, 4.0, 6.0, 3.5], "auto", 3.0) == [0, 0, 1, 1, 1]
    assert scale_groups(values, "group", labels=["a", "a", "b", "b", "b"]) == [0, 0, 1, 1, 1]
    assert scale_groups([5.0, 0.0, float("nan")], "auto") == [0, 1, 2]
    with pytest.raises(ValueError):
        scale_groups(values, "bogus")


# --- colour (z) scale -----------------------------------------------------------------------

def test_zscale_auto_groups_comparable_tracks() -> None:
    m = _matrix([16, 18, 4, 5])
    lims = heatmap_limits(_blocks(m), HeatmapPlotStyle(), m.meta_groups)
    assert lims[0] == pytest.approx(lims[1])
    assert lims[2] == pytest.approx(lims[3])
    assert lims[2][1] < lims[0][1] / 2


def test_zscale_shared_independent_group() -> None:
    m = _matrix([16, 18, 4, 5], groups=["chip", "chip", "cr", "cr"])
    blocks = _blocks(m)
    shared = heatmap_limits(blocks, HeatmapPlotStyle(zscale="shared"), m.meta_groups)
    assert len({lim for lim in shared}) == 1
    indep = heatmap_limits(blocks, HeatmapPlotStyle(zscale="independent"), m.meta_groups)
    assert len({lim for lim in indep}) == 4
    grouped = heatmap_limits(blocks, HeatmapPlotStyle(zscale="group"), m.meta_groups)
    assert grouped[0] == grouped[1] and grouped[2] == grouped[3] and grouped[0] != grouped[2]


def test_single_colorbar_forces_shared() -> None:
    m = _matrix([16, 4])
    lims = heatmap_limits(_blocks(m), HeatmapPlotStyle(independent_scale=False), m.meta_groups)
    assert lims[0] == lims[1]


def test_zmin_zmax_scalar_and_per_sample() -> None:
    m = _matrix([16, 4, 5])
    blocks = _blocks(m)
    fixed = heatmap_limits(blocks, HeatmapPlotStyle(zmin=0, zmax=9, autoscale=False))
    assert fixed == [(0.0, 9.0)] * 3
    per = heatmap_limits(
        blocks, HeatmapPlotStyle(zmin=[0, 0, 0], zmax=[20, None, 6], autoscale=False)
    )
    assert per[0] == (0.0, 20.0) and per[2] == (0.0, 6.0)
    assert 1 < per[1][1] < 20  # automatic (grouped) maximum


# --- profile y axis -------------------------------------------------------------------------

def test_profile_auto_shares_comparable_panels(monkeypatch, tmp_path: Path) -> None:
    got = _draw(monkeypatch, _matrix([8, 7]), HeatmapPlotStyle(), tmp_path)
    assert got["limits"][0] == pytest.approx(got["limits"][1])
    assert got["labels"] == [True, False]


def test_profile_auto_groups_and_labels(monkeypatch, tmp_path: Path) -> None:
    got = _draw(monkeypatch, _matrix([8, 7, 1.5, 1.4]), HeatmapPlotStyle(), tmp_path)
    lim = got["limits"]
    assert lim[0] == pytest.approx(lim[1]) and lim[2] == pytest.approx(lim[3])
    assert lim[2][1] < lim[0][1] / 2
    assert got["labels"] == [True, False, True, False]


def test_profile_forced_shared(monkeypatch, tmp_path: Path) -> None:
    got = _draw(monkeypatch, _matrix([8, 1.5]), HeatmapPlotStyle(profile_yaxis="shared"), tmp_path)
    assert got["limits"][0] == pytest.approx(got["limits"][1])
    assert got["labels"] == [True, False]


def test_profile_user_limits(monkeypatch, tmp_path: Path) -> None:
    got = _draw(
        monkeypatch, _matrix([8, 1.5]), HeatmapPlotStyle(profile_ymin=0, profile_ymax=10), tmp_path
    )
    assert got["limits"] == [pytest.approx((0, 10))] * 2
    assert got["labels"] == [True, False]
    got = _draw(
        monkeypatch, _matrix([8, 1.5]),
        HeatmapPlotStyle(profile_ymin=[0, 0], profile_ymax=[None, 3]), tmp_path,
    )
    assert got["limits"][0][0] == 0 and got["limits"][1] == pytest.approx((0, 3))
    assert got["labels"] == [True, True]


def test_cli_float_list() -> None:
    assert _float_list("") is None
    assert _float_list("2.5") == [2.5]
    assert _float_list("0,auto,3") == [0.0, None, 3.0]


# --- locus tracks ---------------------------------------------------------------------------

def test_track_auto_scale_groups_comparable_tracks() -> None:
    from types import SimpleNamespace

    from genomecanvas.render.renderer import _limits

    tracks = [SimpleNamespace(id=k, type="signal", style={}, scale_group=None) for k in "abc"]
    values = {"a": np.array([0.0, 30.0]), "b": np.array([0.0, 20.0]), "c": np.array([0.0, 1.0])}
    lims = _limits(SimpleNamespace(scales={"mode": "auto"}), tracks, values)
    assert lims["a"] == lims["b"] and lims["a"][1] == 30.0
    assert lims["c"][1] == 1.0
    shared = _limits(SimpleNamespace(scales={"mode": "auto_highest"}), tracks, values)
    assert shared["c"][1] == 30.0
