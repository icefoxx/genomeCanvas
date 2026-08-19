"""Render a real-data GenomeCanvas comparison against the public demo inputs."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from genomecanvas.io.bigwig import read_bigwig
from genomecanvas.io.cytoband import Cytoband, read_cytobands
from genomecanvas.render.plotting import make_grid, plot_signal, save_figure, style_axis

ROOT = Path(__file__).parents[2]
DATA = ROOT / "tests/fixtures/reference_demo"
OUT = Path(__file__).parent
LOCI = [
    ("chr7", 11_002_745, 11_013_547),
    ("chr7", 11_158_374, 11_169_159),
    ("chr7", 11_303_378, 11_313_682),
    ("chr7", 11_394_915, 11_405_318),
]
SIGNALS = [
    ("Input rep1", "GSM5746910_MS_2cell_Input_rep1.bigWig", "#d0d0d0"),
    ("Input rep2", "GSM5746911_MS_2cell_Input_rep2.bigWig", "#d0d0d0"),
    ("IP rep1", "GSM5746912_MS_2cell_IP_rep1.bigWig", "#606060"),
    ("IP rep2", "GSM5746913_MS_2cell_IP_rep2.bigWig", "#606060"),
]
Y_SCALE_MODE = "shared_row"
Y_LIMITS: dict[str, tuple[float, float] | None] = {label: None for label, _, _ in SIGNALS}
SHOW_GENE_TRACK = True
SHOW_CYTOBAND = True
X_LABEL = "GRCm38 chromosome 7 coordinate"
FIGURE_WIDTH = 13.0
PANEL_HEIGHT = 10.0
TRACK_GAP = 0.015
COLUMN_GAP = 0.12
GENE_TRACK_HEIGHT = 3.2
CYTOBAND_HEIGHT = 0.45
GENE_TRACK_SCALE = 1.0
CYTOBAND_SCALE = 1.0
GENE_CYTOBAND_GAP = 0.02
CYTOBAND_MARKER_COLOR = "#e33d32"
CYTOBAND_MARKER_SHAPE = "bar"
GENE_MODEL_LABEL = "none"
TRANSCRIPT_LABEL = "id"


def query(path: Path, chrom: str, start: int, end: int) -> np.ndarray:
    """Query one BigWig interval without loading whole chromosomes."""
    values = np.asarray(read_bigwig(str(path), chrom, start, end), dtype=float)
    return np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)


def genes_for_loci() -> list[list[tuple[int, int, str, str | None]]]:
    """Extract transcript/exon spans overlapping the four displayed loci."""
    results: list[list[tuple[int, int, str, str | None]]] = [[] for _ in LOCI]
    pattern = re.compile(r'transcript_id "([^"]+)"')
    with gzip.open(DATA / "Mus_musculus.GRCm38.102.gtf.gz", "rt", encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"exon", "CDS"}:
                continue
            if fields[0] not in {"7", "chr7"}:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            match = pattern.search(fields[8])
            if match is None:
                continue
            transcript = match.group(1)
            strand = fields[6] if fields[6] in {"+", "-"} else None
            for index, (_, view_start, view_end) in enumerate(LOCI):
                if start < view_end and end > view_start:
                    results[index].append((start, end, transcript, strand))
    return results


def draw_gene(
    axis: plt.Axes, records: list[tuple[int, int, str, str | None]], start: int, end: int,
    *, transcript_label: str = "id", gene_label: str | None = None,
    label_position: str = "bottom",
) -> None:
    """Draw grouped transcript structures with introns, exons, strand, and labels."""
    transcripts: dict[str, list[tuple[int, int, str | None]]] = {}
    for record_start, record_end, transcript, strand in records:
        transcripts.setdefault(transcript, []).append((record_start, record_end, strand))
    spans = sorted(
        ((min(item[0] for item in exons), max(item[1] for item in exons), transcript, exons)
         for transcript, exons in transcripts.items()),
        key=lambda item: (item[0], item[1], item[2]),
    )
    row_ends: list[int] = []
    placed: list[tuple[int, int, str, list[tuple[int, int, str | None]], int]] = []
    for transcript_start, transcript_end, transcript, exons in spans:
        # One transcript per row keeps long stable IDs readable and prevents
        # labels from colliding with neighboring structures.
        row = len(row_ends)
        row_ends.append(transcript_end)
        placed.append((transcript_start, transcript_end, transcript, exons, row))
    for transcript_start, transcript_end, transcript, exons, row in placed:
        y = len(row_ends) - row
        left, right = max(transcript_start, start), min(transcript_end, end)
        if left >= right:
            continue
        strand = next((item[2] for item in exons if item[2] in {"+", "-"}), None)
        axis.plot([left, right], [y, y], color="#777777", linewidth=0.65)
        for exon_start, exon_end, _ in exons:
            exon_left, exon_right = max(exon_start, start), min(exon_end, end)
            if exon_left < exon_right:
                axis.add_patch(
                    plt.Rectangle(
                        (exon_left, y - 0.13), exon_right - exon_left, 0.26, color="#222222"
                    )
                )
        direction = 1 if strand != "-" else -1
        strand_color = "#1f78b4" if direction > 0 else "#d73027"
        arrow_x = right - 0.08 * (right - left) if direction > 0 else left + 0.08 * (right - left)
        axis.annotate(
            "", xy=(arrow_x + direction * 0.04 * (right - left), y), xytext=(arrow_x, y),
            arrowprops={
                "arrowstyle": "-|>", "color": strand_color, "lw": 1.4,
                "mutation_scale": 9,
            },
        )
        if transcript_label != "none":
            label_y = y + 0.18 if label_position == "top" else y - 0.34
            label_va = "bottom" if label_position == "top" else "top"
            axis.text(
                (left + right) / 2, label_y, transcript, fontsize=5.5, color="#555555",
                ha="center", va=label_va, clip_on=False,
            )
    axis.set_xlim(start, end)
    axis.set_ylim(-0.55, max(1, len(row_ends)) + 0.3)
    axis.set_yticks([])
    axis.spines[["top", "right", "left", "bottom"]].set_visible(False)
    if gene_label:
        axis.text(0.01, 1.02, gene_label, transform=axis.transAxes, fontsize=7, va="bottom")


def draw_cytoband(
    axis: plt.Axes, bands: list[Cytoband], marker_start: int, marker_end: int,
    marker_style: str = CYTOBAND_MARKER_SHAPE, marker_color: str = CYTOBAND_MARKER_COLOR,
) -> None:
    """Draw actual chromosome bands and map the locus onto the ideogram."""
    chromosome_end = bands[-1].end
    stain_colors = {
        "gneg": "#f7f7f7", "gpos25": "#d9d9d9", "gpos33": "#bdbdbd",
        "gpos50": "#969696", "gpos66": "#737373", "gpos75": "#525252",
        "gpos100": "#111111", "acen": "#e6614f", "gvar": "#d9d9d9", "stalk": "#bdbdbd",
    }
    for band in bands:
        axis.add_patch(
            plt.Rectangle(
                (band.start / chromosome_end, 0.3), (band.end - band.start) / chromosome_end,
                0.4, color=stain_colors.get(band.stain, "#bdbdbd"), ec="#555", lw=0.2,
            )
        )
    marker_position = ((marker_start + marker_end) / 2) / chromosome_end
    if marker_style == "triangle":
        axis.scatter([marker_position], [0.1], marker="v", color=marker_color, s=60)
    elif marker_style == "line":
        axis.plot([marker_position, marker_position], [0.2, 0.8], color=marker_color, linewidth=2.5)
    elif marker_style == "bar":
        axis.add_patch(
            plt.Rectangle((marker_position - 0.008, 0.22), 0.016, 0.56, color=marker_color)
        )
    elif marker_style in {"circle", "diamond"}:
        marker = "o" if marker_style == "circle" else "D"
        axis.scatter([marker_position], [0.5], marker=marker, color=marker_color, s=45)
    else:
        raise ValueError("marker_style must be bar, line, or triangle")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")


def main() -> None:
    genes = genes_for_loci()
    cytobands = read_cytobands(DATA / "mm10_cytoBandIdeo.txt.gz", chrom="chr7")
    shared_limits: dict[str, tuple[float, float]] = {}
    for label, filename, _ in SIGNALS:
        maxima = [
            float(query(DATA / filename, chrom, start, end).max())
            for chrom, start, end in LOCI
        ]
        explicit = Y_LIMITS[label]
        shared_limits[label] = explicit or (0.0, max(maxima) * 1.05)
    figure, axes = make_grid(
        6, 4, figsize=(FIGURE_WIDTH, PANEL_HEIGHT),
        height_ratios=[2, 2, 2, 2, 0.8, 0.35],
    )
    figure.subplots_adjust(hspace=TRACK_GAP, wspace=0.05, top=0.9, bottom=0.06)
    for column, (chrom, start, end) in enumerate(LOCI):
        for row, (label, filename, color) in enumerate(SIGNALS):
            values = query(DATA / filename, chrom, start, end)
            positions = np.linspace(start, end, len(values), endpoint=False)
            axis = axes[row, column]
            plot_signal(axis, positions, values, color=color, alpha=1.0, linewidth=0)
            if Y_SCALE_MODE == "shared_row":
                axis.set_ylim(*shared_limits[label])
            elif Y_LIMITS[label] is not None:
                axis.set_ylim(*Y_LIMITS[label])
            axis.set_xlim(start, end)
            style_axis(axis)
            axis.spines[["top", "right", "left"]].set_color("#d0d0d0")
            axis.tick_params(labelsize=7)
            if column == 0:
                axis.set_ylabel(label, fontsize=8)
            if row == 0:
                axis.set_title(
                    f"real GRCm38 locus | {chrom}:{start:,}-{end:,}",
                    fontsize=8, weight="bold",
                )
        if SHOW_GENE_TRACK:
            draw_gene(axes[4, column], genes[column], start, end)
        else:
            axes[4, column].set_visible(False)
        if SHOW_CYTOBAND:
            draw_cytoband(axes[5, column], cytobands, start, end)
        else:
            axes[5, column].set_visible(False)
        axes[3, column].set_xlabel(X_LABEL)
        for row in range(5):
            axes[row, column].set_xlim(start, end)
    figure.suptitle("GenomeCanvas real-data faceted comparison", fontsize=13, weight="bold")
    save_figure(figure, OUT, "p01_genomecanvas_realdata")


if __name__ == "__main__":
    main()
