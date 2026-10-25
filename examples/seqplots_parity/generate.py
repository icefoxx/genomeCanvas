"""Render GenomeCanvas figures from the official SeqPlots ce10 example files."""

from __future__ import annotations

from pathlib import Path

from genomecanvas.io.bed import read_bed
from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.pipeline.matrix import compute_region_matrix, sort_matrix
from genomecanvas.pipeline.metadata import SampleMeta
from genomecanvas.pipeline.plot_style import HeatmapPlotStyle, ProfilePlotStyle
from genomecanvas.pipeline.plotting import plot_enrichment_heatmap, plot_profile_lines
from genomecanvas.render.renderer import render

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests/fixtures/seqplots_demo"
OUT = Path(__file__).parent / "output"
OUT.mkdir(exist_ok=True)

BW = DATA / "GSM1208360_chrI_100Kb_q5_sample.bw"
PEAKS = DATA / "GSM1208361_chrI_100Kb_PeakCalls.bed"
TX = DATA / "Transcripts_ce10_chrI_100Kb.bed"
GTF = DATA / "transcripts_ce10_chrI_100Kb.gtf"


def bed12_to_gtf(bed: Path, gtf: Path) -> Path:
    """Write a GTF with exon blocks from a BED12 (or BED6) file."""
    lines: list[str] = []
    for feature in read_bed(bed):
        interval = feature.interval
        name = feature.name or "tx"
        strand = interval.strand or "."
        gene = name.rsplit(".", 1)[0]
        attrs = f'gene_id "{gene}"; transcript_id "{name}"; gene_name "{gene}";'
        lines.append(
            f"{interval.chrom}\tseqplots\texon\t{interval.start + 1}\t{interval.end}\t.\t"
            f"{strand}\t.\t{attrs}"
        )
    gtf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return gtf


def sample() -> SampleMeta:
    return SampleMeta(str(BW), "H3K4me3", "chip", "h3k4me3", "Blues", "#1f4e79")


def heatmap(
    beds: list[tuple[str, str]],
    output: Path,
    style: HeatmapPlotStyle,
    *,
    ref_point: str = "center",
) -> None:
    matrix = compute_region_matrix(
        [sample()], beds, flank=500, binsize=25, ref_point=ref_point,  # type: ignore[arg-type]
        skip_zeros=True, engine="auto",
    )
    matrix = sort_matrix(matrix, method="mean")
    plot_enrichment_heatmap(matrix, output, style=style)
    plot_profile_lines(
        matrix,
        output.with_name(output.stem.replace("heatmap", "profile") + output.suffix),
        style=ProfilePlotStyle(
            title=style.title,
            labels=style.labels,
            xlab=style.xlab,
            ylab="mean coverage",
        ),
    )


def main() -> None:
    bed12_to_gtf(TX, GTF)
    spec = FigureSpec.parse_obj(
        {
            "genome": {"assembly": "ce10"},
            "region": {"chrom": "chrI", "start": 0, "end": 100_000},
            "gene_model": str(GTF),
            "collapse_transcripts": True,
            "tracks": [
                {
                    "id": "h3k4me3",
                    "type": "signal",
                    "label": "H3K4me3",
                    "source": {"path": str(BW)},
                    "style": {"color": "#1f4e79", "ymin": 0},
                },
                {
                    "id": "peaks",
                    "type": "peak",
                    "label": "peaks",
                    "source": {"path": str(PEAKS)},
                    "style": {"color": "#c0392b", "show_labels": False},
                },
                {"id": "genes", "type": "gene"},
            ],
            "layout": {
                "width": 10,
                "height": 5.2,
                "x_tick_unit": "kb",
                "transcript_label": "name",
            },
            "style": {"title": "SeqPlots example · C. elegans chrI 100 kb (ce10)"},
        }
    )
    render(spec, OUT / "seqplots_locus.png")

    heatmap(
        [(str(TX), "Transcripts")],
        OUT / "seqplots_transcript_tss_heatmap.png",
        HeatmapPlotStyle(
            title="H3K4me3 at transcript TSS (SeqPlots-style anchored)",
            labels=("H3K4me3",),
            xlab="distance from TSS",
            ylab="transcripts",
            colorbar_label="coverage",
            height=6.5,
        ),
        ref_point="TSS",
    )
    heatmap(
        [(str(TX), "Transcripts")],
        OUT / "seqplots_transcript_heatmap.png",
        HeatmapPlotStyle(
            title="H3K4me3 at transcripts (SeqPlots example)",
            labels=("H3K4me3",),
            xlab="distance from transcript center",
            ylab="transcripts",
            colorbar_label="coverage",
            height=6.5,
        ),
    )
    heatmap(
        [(str(PEAKS), "Peaks")],
        OUT / "seqplots_peak_heatmap.png",
        HeatmapPlotStyle(
            title="H3K4me3 at peak calls (SeqPlots example)",
            labels=("H3K4me3",),
            xlab="distance from peak center",
            ylab="peaks",
            colorbar_label="coverage",
            height=6.0,
        ),
    )
    heatmap(
        [(str(TX), "Transcripts"), (str(PEAKS), "Peaks")],
        OUT / "seqplots_features_heatmap.png",
        HeatmapPlotStyle(
            title="H3K4me3 · transcripts vs peaks (SeqPlots getPlotSetArray)",
            labels=("H3K4me3",),
            xlab="distance from feature center",
            ylab="features",
            colorbar_label="coverage",
            height=7.0,
        ),
    )
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
