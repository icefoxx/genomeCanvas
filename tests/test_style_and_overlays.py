from pathlib import Path

from genomecanvas.interact import FigureEdit, OverlayAnnotation, apply_edit, hit_at_genome
from genomecanvas.interact.hits import HitRegistry
from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.model.overlays import HitTarget
from genomecanvas.render.renderer import render

FIXTURE = Path(__file__).parent / "fixtures/mini_genome"


def _spec(**extra):
    payload = {
        "genome": {"assembly": "synthetic"},
        "region": {"chrom": "chr1", "start": 0, "end": 300},
        "gene_model": str(FIXTURE / "mini.gtf"),
        "style": {
            "title": "Styled figure",
            "font_family": "DejaVu Sans",
            "track_label_font": {"size": 9, "style": "italic", "color": "#113355"},
            "gene_label_font": {"size": 8, "style": "italic"},
        },
        "tracks": [
            {
                "id": "signal", "type": "signal", "label": "A",
                "source": {"path": str(FIXTURE / "signal_a.bedgraph")},
                "height": 1.4,
                "style": {"color": "#d35400", "font_style": "italic"},
            },
            {
                "id": "genes", "type": "gene", "height": 0.8,
                "style": {"label_mode": "name"},
            },
            {
                "id": "iso", "type": "isoform", "height": 1.0,
                "style": {"collapse_transcripts": False},
            },
            {
                "id": "peaks", "type": "peak",
                "source": {"path": str(FIXTURE / "peaks.bed")},
            },
            {
                "id": "motifs", "type": "motif",
                "source": {"path": str(FIXTURE / "motifs.bed")},
            },
            {
                "id": "loops", "type": "loop",
                "source": {"path": str(FIXTURE / "loops.bedpe")},
            },
            {
                "id": "snvs", "type": "variant",
                "source": {"path": str(FIXTURE / "variants.vcf")},
            },
        ],
        "overlays": [
            {
                "id": "note1", "kind": "label", "start": 50, "end": 80,
                "text": "peak", "track_id": "signal", "font_style": "italic",
            }
        ],
    }
    payload.update(extra)
    return FigureSpec.parse_obj(payload)


def test_render_all_track_types_and_hits(tmp_path: Path) -> None:
    output = tmp_path / "all.png"
    render(_spec(), output)
    assert output.is_file()
    hits_path = tmp_path / "all.hits.json"
    assert hits_path.is_file()
    text = hits_path.read_text(encoding="utf-8")
    assert "track:signal" in text
    assert "overlay:note1" in text
    assert "gene:" in text


def test_apply_edit_adds_and_removes_overlay() -> None:
    spec = _spec()
    added = apply_edit(
        spec,
        FigureEdit(
            action="add_overlay",
            overlay=OverlayAnnotation(id="n2", kind="text", start=100, text="hello"),
        ),
    )
    assert overlay_ids(added) == ["note1", "n2"]
    removed = apply_edit(added, FigureEdit(action="remove_overlay", overlay_id="note1"))
    assert overlay_ids(removed) == ["n2"]
    taller = apply_edit(
        spec, FigureEdit(action="set_track_height", track_id="signal", height=2.5)
    )
    assert taller.tracks[0].height == 2.5
    italic = apply_edit(
        spec,
        FigureEdit(
            action="update_track_style",
            track_id="signal",
            style={"font_style": "italic", "font_size": 12},
        ),
    )
    assert italic.tracks[0].style["font_size"] == 12
    titled = apply_edit(spec, FigureEdit(action="set_title", title="Clicked title"))
    assert titled.style.title == "Clicked title"


def test_hit_at_genome_prefers_later_overlay() -> None:
    registry = HitRegistry()
    registry.add(
        HitTarget(
            id="t", kind="track", track_id="signal", chrom="chr1",
            start=0, end=300, axes_index=0, x=0, y=0,
        )
    )
    registry.add(
        HitTarget(
            id="o", kind="overlay", overlay_id="note1", track_id="signal",
            chrom="chr1", start=50, end=80, axes_index=0, x=50, y=1, text="peak",
        )
    )
    hit = hit_at_genome(registry.targets, "chr1", 60, track_id="signal")
    assert hit is not None
    assert hit.kind == "overlay"


def overlay_ids(spec: FigureSpec) -> list[str]:
    return [item.id for item in spec.overlays]
