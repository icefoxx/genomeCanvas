from pathlib import Path

from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.render.renderer import render


def test_render_writes_required_artifacts(tmp_path: Path) -> None:
    fixture = Path(__file__).parent / "fixtures/mini_genome/signal_a.bedgraph"
    spec = FigureSpec.parse_obj({
        "genome": {"assembly": "synthetic"},
        "region": {"chrom": "chr1", "start": 0, "end": 300},
        "tracks": [{"id": "signal", "type": "signal", "source": {"path": str(fixture)}}],
    })
    output = tmp_path / "figure.pdf"
    render(spec, output)
    assert output.is_file()
    assert (tmp_path / "figure.svg").is_file()
    assert (tmp_path / "figure.spec.yaml").is_file()
    assert (tmp_path / "figure.manifest.json").is_file()
