"""Smoke tests for TASK-001 repository infrastructure."""

import genomecanvas
from genomecanvas.core.coordinates import internal_to_one_based_closed, one_based_closed_to_internal
from genomecanvas.core.intervals import GenomicInterval
from genomecanvas.layout.interval_packing import pack_intervals
from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.model.region import GenomeRegion


def test_package_imports() -> None:
    assert genomecanvas.__version__ == "0.1.0"


def test_interval_operations() -> None:
    interval = GenomicInterval("chr1", 10, 20)
    assert interval.length == 10
    assert interval.contains(10)
    assert not interval.contains(20)
    assert interval.intersection(GenomicInterval("chr1", 15, 25)) == GenomicInterval("chr1", 15, 20)
    assert interval.distance_to(GenomicInterval("chr1", 25, 30)) == 5


def test_coordinate_round_trip() -> None:
    internal = one_based_closed_to_internal(11, 20)
    assert internal == (10, 20)
    assert internal_to_one_based_closed(*internal) == (11, 20)


def test_region_parser() -> None:
    region = GenomeRegion.parse("chr7:55,000,000-55,100,000", assembly="hg38")
    assert region.start == 55_000_000
    assert region.end == 55_100_000


def test_interval_packing_is_stable_and_non_overlapping() -> None:
    intervals = ((0, 10), (5, 8), (10, 15), (2, 4))
    assert pack_intervals(intervals) == ((0, 2), (3, 1))


def test_figure_spec_round_trip() -> None:
    spec = FigureSpec.parse_obj({
        "genome": {"assembly": "hg38"},
        "region": {"chrom": "chr7", "start": 10, "end": 20},
        "tracks": [{"id": "signal", "type": "signal", "source": {"path": "x.bg"}}],
    })
    loaded = FigureSpec.parse_raw(spec.json())
    assert loaded == spec
