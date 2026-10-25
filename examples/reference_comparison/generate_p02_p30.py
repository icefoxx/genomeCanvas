"""Generate real-data p02-p30 behavior variants independently.

The filenames mirror the public demo's behavior categories only. Rendering is
implemented here with GenomeCanvas's BigWig adapter, explicit scale policies,
and the downloaded GRCm38/mm10 data; no reference plotting source is imported.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
from matplotlib.ticker import FuncFormatter, LinearLocator, NullFormatter, NullLocator

from genomecanvas.io.cytoband import read_cytobands
from genomecanvas.render.plotting import make_grid, plot_signal, save_figure, style_axis

from .generate import (
    COLUMN_GAP,
    CYTOBAND_HEIGHT,
    CYTOBAND_MARKER_COLOR,
    CYTOBAND_MARKER_SHAPE,
    CYTOBAND_SCALE,
    DATA,
    FIGURE_WIDTH,
    GENE_CYTOBAND_GAP,
    GENE_MODEL_LABEL,
    GENE_TRACK_HEIGHT,
    GENE_TRACK_SCALE,
    PANEL_HEIGHT,
    SIGNALS,
    TRACK_GAP,
    TRANSCRIPT_LABEL,
    draw_cytoband,
    draw_gene,
    query,
)

OUT = Path(__file__).parent / "output"
P03_LOCUS = ("chr14", 54_832_580, 55_004_465)
P03_TRACK_COLORS = ["#b3cde3", "#6497b1", "#005b96", "#03396c"]
P22_LOCI = [
    ("chr7", 10_897_740, 10_908_050),
    ("chr7", 11_002_745, 11_013_547),
]
P22_HIGHLIGHTS = [
    (10_904_000, 10_906_000, "#e88989"),
    (11_008_000, 11_010_000, "#8f8fe8"),
]
REGION_ZOOM_LOCI = [
    ("chr14", 54_835_580, 55_001_465),
    ("chr19", 5_790_000, 5_810_000),
]
LOCI7 = [
    ("chr7", 11_002_745, 11_013_547),
    ("chr7", 11_158_374, 11_169_159),
    ("chr7", 11_303_378, 11_313_682),
    ("chr7", 11_394_915, 11_405_318),
]
LOCI14 = [
    ("chr14", 54_938_921, 54_969_927),
    ("chr14", 54_967_684, 54_997_626),
    ("chr14", 54_880_377, 54_891_234),
    ("chr14", 54_889_500, 54_901_169),
]


def records_for_loci(
    loci: list[tuple[str, int, int]],
    gene_names: list[str] | None = None,
) -> list[list[tuple[int, int, str, str | None]]]:
    """Extract real GTF exon/CDS records for arbitrary displayed loci."""
    result: list[list[tuple[int, int, str, str | None]]] = [[] for _ in loci]
    chromosomes = {chrom.removeprefix("chr") for chrom, _, _ in loci}
    pattern = re.compile(r'transcript_id "([^"]+)"')
    gene_pattern = re.compile(r'gene_name "([^"]+)"')
    with gzip.open(DATA / "Mus_musculus.GRCm38.102.gtf.gz", "rt", encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0].removeprefix("chr") not in chromosomes:
                continue
            if fields[2] not in {"exon", "CDS"}:
                continue
            match = pattern.search(fields[8])
            if match is None:
                continue
            gene_match = gene_pattern.search(fields[8])
            start, end = int(fields[3]) - 1, int(fields[4])
            strand = fields[6] if fields[6] in {"+", "-"} else None
            for index, (chrom, view_start, view_end) in enumerate(loci):
                same_chrom = fields[0].removeprefix("chr") == chrom.removeprefix("chr")
                selected_gene = gene_names is None or (
                    gene_match is not None and gene_match.group(1) == gene_names[index]
                )
                if same_chrom and selected_gene and start < view_end and end > view_start:
                    result[index].append((start, end, match.group(1), strand))
    return result


def panel_grid(
    name: str,
    loci: list[tuple[str, int, int]],
    *,
    gene: bool = True,
    cyto: bool = True,
    compact: bool = False,
    box: bool = False,
    colors: list[str] | None = None,
    fixed: tuple[float, float] | None = None,
    gene_top: bool = False,
    scale_mode: str = "auto",
    gene_names: list[str] | None = None,
    aspect_ratio: float | None = None,
    highlights: list[tuple[int, int, str]] | None = None,
    gene_track_scale: float = GENE_TRACK_SCALE,
    cytoband_scale: float = CYTOBAND_SCALE,
    gene_cytoband_gap: float = GENE_CYTOBAND_GAP,
    cytoband_marker_shape: str = CYTOBAND_MARKER_SHAPE,
    cytoband_marker_color: str = CYTOBAND_MARKER_COLOR,
    gene_model_label: str = GENE_MODEL_LABEL,
    transcript_label: str = TRANSCRIPT_LABEL,
    highlight_gene: bool = False,
    upper_annotation_ratio: float | None = None,
    y_axis_mode: str = "all",
    y_value_transform: str = "identity",
    y_divisor: float = 1.0,
    y_tick_count: int = 5,
    y_tick_decimals: int = 1,
    y_tick_label_pad: float = 3.0,
    show_y_ticks: bool = True,
    y_tick_length: float = 3.0,
    show_y_axis: bool = True,
    show_major_frame: bool = True,
) -> None:
    """Render four real-data facets with configurable presentation layers."""
    records = records_for_loci(loci, gene_names=gene_names)
    bands = read_cytobands(DATA / "mm10_cytoBandIdeo.txt.gz", chrom=loci[0][0])
    rows = 6 if gene and cyto else 5 if gene or cyto else 4
    signal_offset = 1 if gene and gene_top else 0
    ncols = len(loci)
    height_ratios = [2] * 4
    if gene:
        height_ratios.insert(0 if gene_top else 4, GENE_TRACK_HEIGHT * gene_track_scale)
    if cyto:
        height_ratios.append(CYTOBAND_HEIGHT * cytoband_scale)
    if upper_annotation_ratio is not None and gene and cyto:
        upper_total = sum(height_ratios[:4])
        lower_total = sum(height_ratios[4:])
        lower_factor = upper_total / (upper_annotation_ratio * lower_total)
        height_ratios[4:] = [height * lower_factor for height in height_ratios[4:]]
    figure, axes = make_grid(
        rows, ncols,
        figsize=(
            FIGURE_WIDTH * ncols / 4,
            FIGURE_WIDTH * aspect_ratio
            if aspect_ratio is not None
            else PANEL_HEIGHT * (0.8 if compact else 1.0),
        ),
        constrained_layout=False,
        gridspec_kw={"height_ratios": height_ratios},
    )
    figure.subplots_adjust(
        hspace=max(TRACK_GAP, gene_cytoband_gap), wspace=COLUMN_GAP,
        top=0.9, bottom=0.06,
    )
    colors = colors or [item[2] for item in SIGNALS]
    matrix = np.array([
        [
            float(
                _display_values(
                    query(DATA / SIGNALS[row][1], *locus), y_value_transform, y_divisor
                ).max()
            )
            for row in range(4)
        ]
        for locus in loci
    ])
    if scale_mode == "auto":
        column_waste = sum(float(column.max() * 4 - column.sum()) for column in matrix)
        row_waste = sum(float(row.max() * len(loci) - row.sum()) for row in matrix.T)
        scale_mode = "column" if column_waste <= row_waste else "row"
    for column, (chrom, start, end) in enumerate(loci):
        for row, (label, filename, _) in enumerate(SIGNALS):
            values = _display_values(
                query(DATA / filename, chrom, start, end), y_value_transform, y_divisor
            )
            positions = np.linspace(start, end, len(values), endpoint=False)
            axis = axes[signal_offset + row, column]
            plot_signal(axis, positions, values, color=colors[row], alpha=1.0, linewidth=0)
            if highlights is not None:
                highlight_start, highlight_end, highlight_color = highlights[column]
                axis.axvspan(
                    highlight_start, highlight_end, color=highlight_color,
                    alpha=0.32, zorder=0,
                )
            axis.set_xlim(start, end)
            if fixed is not None:
                axis.set_ylim(*fixed)
            elif scale_mode == "column":
                axis.set_ylim(0, float(matrix[column].max()) * 1.05)
            else:
                axis.set_ylim(0, float(matrix[:, row].max()) * 1.05)
            style_axis(axis)
            axis.spines[["top", "right", "left"]].set_color("#d0d0d0")
            axis.tick_params(labelsize=7)
            show_left = y_axis_mode == "all" or y_axis_mode == "left" or (
                y_axis_mode == "outer" and column == 0
            )
            show_right = y_axis_mode == "right" or (
                y_axis_mode == "outer" and column == len(loci) - 1
            )
            if show_left:
                if column == 0:
                    axis.set_ylabel(label, fontsize=8)
                axis.yaxis.set_ticks_position("left")
                axis.tick_params(axis="y", left=True, labelleft=True)
            elif show_right:
                axis.yaxis.set_ticks_position("right")
                axis.tick_params(axis="y", left=False, labelleft=False, right=True, labelright=True)
            else:
                axis.set_yticks([])
                axis.tick_params(
                    axis="y", left=False, right=False, labelleft=False, labelright=False
                )
            axis.spines["left"].set_visible(column == 0)
            axis.spines["right"].set_visible(column == len(loci) - 1)
            if not show_major_frame:
                for spine in axis.spines.values():
                    spine.set_visible(False)
            axis.yaxis.set_major_locator(LinearLocator(y_tick_count))
            axis.yaxis.set_major_formatter(
                FuncFormatter(lambda value, _: f"{value:.{y_tick_decimals}f}")
            )
            axis.tick_params(
                axis="y", length=y_tick_length if show_y_ticks else 0,
                pad=y_tick_label_pad,
            )
            if not show_y_axis:
                axis.tick_params(
                    axis="y", left=False, right=False, labelleft=False, labelright=False
                )
                axis.yaxis.set_visible(False)
            if row == 0:
                axis.set_title(f"{chrom}:{start:,}-{end:,}", fontsize=8, weight="bold")
            axis.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
            axis.set_xticks([])
            axis.xaxis.set_major_locator(NullLocator())
            axis.xaxis.set_minor_locator(NullLocator())
            axis.xaxis.set_major_formatter(NullFormatter())
            axis.xaxis.set_minor_formatter(NullFormatter())
            axis.xaxis.set_visible(False)
            axis.set_xticklabels([])
            if box:
                for spine in axis.spines.values():
                    spine.set_visible(True)
                    spine.set_color("#bdbdbd")
        if gene:
            target = 0 if gene_top else 4
            gene_name = gene_names[column] if gene_names else f"locus {column + 1}"
            draw_gene(
                axes[target, column], records[column], start, end,
                transcript_label=transcript_label,
                gene_label=gene_name if gene_model_label != "none" else None,
            )
            axes[target, column].xaxis.set_visible(False)
            if highlights is not None and highlight_gene:
                highlight_start, highlight_end, highlight_color = highlights[column]
                axes[target, column].axvspan(
                    highlight_start, highlight_end, color=highlight_color,
                    alpha=0.32, zorder=0,
                )
            isoforms = len({record[2] for record in records[column]})
            axes[0, column].text(
                0.5, 1.12, f"{gene_name} | {loci[column][0]} | {isoforms} isoforms",
                transform=axes[0, column].transAxes, ha="center", va="bottom",
                fontsize=9, fontstyle="italic",
            )
        if cyto:
            target = 5 if gene else 4
            draw_cytoband(
                axes[target, column], bands, start, end,
                marker_style=cytoband_marker_shape,
                marker_color=cytoband_marker_color,
            )
            axes[target, column].xaxis.set_visible(False)
        if not gene:
            axes[3, column].set_xlabel(f"{chrom}:{start:,}-{end:,}", fontsize=7)
    figure.suptitle(f"GenomeCanvas real-data panel {name}", fontsize=12, weight="bold")
    save_figure(figure, OUT, name, dpi=150)


def _display_values(values: np.ndarray, transform: str, divisor: float) -> np.ndarray:
    """Apply a display-only y-axis transform."""
    if transform == "divide":
        return values / divisor
    if transform == "log2":
        return np.log2(np.maximum(values, np.finfo(float).tiny))
    if transform == "log10":
        return np.log10(np.maximum(values, np.finfo(float).tiny))
    return values


def region_panel(
    name: str,
    regions: list[tuple[int, int, str, str]],
    track_colors: list[str] | None = None,
) -> None:
    """Render one wide real-data locus with one or more colored regions."""
    chrom, start, end = P03_LOCUS
    records = records_for_loci([P03_LOCUS])[0]
    bands = read_cytobands(DATA / "mm10_cytoBandIdeo.txt.gz", chrom=chrom)
    figure, axes = make_grid(
        7, 1, figsize=(14, 14), sharex=False, constrained_layout=False,
        gridspec_kw={"height_ratios": [2, 2, 2, 2, 0.65, 6.5, 0.45]},
    )
    axes = axes[:, 0]
    figure.subplots_adjust(hspace=TRACK_GAP, top=0.91, bottom=0.07, left=0.07, right=0.96)
    track_colors = track_colors or P03_TRACK_COLORS
    for row, (label, filename, _) in enumerate(SIGNALS):
        values = query(DATA / filename, chrom, start, end)
        positions = np.linspace(start, end, len(values), endpoint=False)
        axis = axes[row]
        plot_signal(axis, positions, values, color=track_colors[row], alpha=1.0, linewidth=0)
        axis.set_xlim(start, end)
        axis.set_ylim(0, max(1.0, float(values.max()) * 1.05))
        axis.set_ylabel(label, fontsize=9)
        axis.tick_params(axis="x", bottom=False, labelbottom=False)
        style_axis(axis)
        axis.spines[["top", "right", "left"]].set_color("#d0d0d0")
    region_axis = axes[4]
    for region_start, region_end, label, color in regions:
        region_axis.axvspan(region_start, region_end, color=color, alpha=0.42)
        region_axis.text(
            (region_start + region_end) / 2, 0.62, label,
            ha="center", va="center", fontsize=8, clip_on=False,
        )
    region_axis.set_xlim(start, end)
    region_axis.set_ylim(0, 1)
    region_axis.set_yticks([])
    region_axis.spines[["top", "right", "left", "bottom"]].set_visible(False)
    draw_gene(axes[5], records, start, end)
    axes[5].xaxis.set_visible(False)
    axes[5].set_xlabel(f"{chrom}:{start:,}-{end:,}", fontsize=8)
    draw_cytoband(axes[6], bands, start, end)
    axes[6].xaxis.set_visible(False)
    figure.suptitle(
        f"{chrom}:{start:,}-{end:,} | GenomeCanvas real-data regions",
        fontsize=13, weight="bold",
    )
    figure.text(
        0.5, 0.925, f"{regions[0][2]} | {chrom} | {len(regions)} region(s)",
        ha="center", fontsize=10, style="italic",
    )
    save_figure(figure, OUT, name, dpi=150)


def trap_panel(
    name: str,
    *,
    gene_top: bool,
    trap_color: str = "#9ecae1",
) -> None:
    """Render one locus with an explicit trap between content and cytoband."""
    chrom, start, end = P22_LOCI[0]
    records = records_for_loci([P22_LOCI[0]])[0]
    bands = read_cytobands(DATA / "mm10_cytoBandIdeo.txt.gz", chrom=chrom)
    rows = 7 if gene_top else 6
    signal_offset = 1 if gene_top else 0
    gene_row = 0 if gene_top else None
    trap_row = 5 if gene_top else 4
    cyto_row = 6 if gene_top else 5
    figure, axes = make_grid(
        rows, 1, figsize=(10, 12), constrained_layout=False,
        gridspec_kw={"height_ratios": ([3.0] if gene_top else []) + [2, 2, 2, 2, 0.45, 0.45]},
    )
    figure.subplots_adjust(hspace=0.015, top=0.91, bottom=0.06, left=0.08, right=0.96)
    signal_arrays = [query(DATA / filename, chrom, start, end) for _, filename, _ in SIGNALS]
    shared_max = max(float(values.max()) for values in signal_arrays) * 1.05
    for row, ((label, _, color), values) in enumerate(zip(SIGNALS, signal_arrays, strict=True)):
        axis = axes[signal_offset + row, 0]
        positions = np.linspace(start, end, len(values), endpoint=False)
        plot_signal(axis, positions, values, color=color, alpha=1.0, linewidth=0)
        axis.set_xlim(start, end)
        axis.set_ylim(0, shared_max)
        axis.set_ylabel(label, fontsize=9)
        axis.tick_params(axis="x", bottom=False, labelbottom=False)
        style_axis(axis)
        axis.spines[["top", "right", "left"]].set_color("#d0d0d0")
    if gene_top and gene_row is not None:
        draw_gene(axes[gene_row, 0], records, start, end)
        axes[gene_row, 0].xaxis.set_visible(False)
        axes[gene_row, 0].set_title(
            f"Zscan4b | {chrom}:{start:,}-{end:,}", fontsize=10, style="italic"
        )
    trap = axes[trap_row, 0]
    chromosome_end = bands[-1].end
    marker_position = ((start + end) / 2) / chromosome_end
    marker_half_width = max(0.018, (end - start) / chromosome_end / 2)
    trap.fill(
        [0.0, 1.0, marker_position + marker_half_width, marker_position - marker_half_width],
        [0.9, 0.9, 0.05, 0.05],
        color=trap_color, alpha=0.85,
    )
    trap.set_xlim(0, 1)
    trap.set_ylim(0, 1)
    trap.axis("off")
    draw_cytoband(axes[cyto_row, 0], bands, start, end)
    figure.suptitle(
        f"{chrom}:{start:,}-{end:,} | GenomeCanvas trap composition",
        fontsize=12, weight="bold",
    )
    save_figure(figure, OUT, name, dpi=150)


def ratio_panel(name: str, loci: list[tuple[str, int, int]], pseudocount: float) -> None:
    """Render explicit log2 ratios for Input rep2 versus Input rep1."""
    figure, axes = make_grid(4, 4, figsize=(13, 8))
    for column, (chrom, start, end) in enumerate(loci):
        left = query(DATA / SIGNALS[0][1], chrom, start, end)
        right = query(DATA / SIGNALS[1][1], chrom, start, end)
        ratio = np.log2((right + pseudocount) / (left + pseudocount))
        positions = np.linspace(start, end, len(ratio), endpoint=False)
        for row in range(4):
            axis = axes[row, column]
            axis.plot(positions, ratio, color="#555555", linewidth=0.45)
            axis.axhline(0, color="#222222", linewidth=0.6)
            axis.set_xlim(start, end)
            axis.set_ylabel("log2 ratio" if column == 0 else "")
            axis.set_title(f"{chrom}:{start:,}-{end:,}" if row == 0 else "", fontsize=8)
    figure.suptitle(
        f"GenomeCanvas log2 ratio · pseudocount {pseudocount}", fontsize=12, weight="bold"
    )
    save_figure(figure, OUT, name, dpi=150)


def zoom_panel(
    name: str, loci: list[tuple[str, int, int]], trap: bool = True,
    zoom_side: str = "bottom", highlight_color: str = "#f2c14e",
) -> None:
    """Render compact parent/gene/trap/child/gene/cytoband compositions."""
    displayed = loci[:2]
    ncols = len(displayed)
    if zoom_side not in {"top", "bottom"}:
        raise ValueError("zoom_side must be 'top' or 'bottom'")
    parent_row_start = 0 if zoom_side == "bottom" else 6
    child_row_start = 6 if zoom_side == "bottom" else 1
    parent_gene_row = 4 if zoom_side == "bottom" else 10
    child_gene_row = 10 if zoom_side == "bottom" else 0
    height_ratios = (
        [2, 2, 2, 2, 0.7, 0.5, 2, 2, 2, 2, 0.7, 0.45]
        if zoom_side == "bottom"
        else [0.7, 2, 2, 2, 2, 0.5, 2, 2, 2, 2, 0.7, 0.45]
    )
    figure, axes = make_grid(
        12, ncols, figsize=(10 * ncols, 14), constrained_layout=False,
        gridspec_kw={"height_ratios": height_ratios},
    )
    figure.subplots_adjust(hspace=0.01, wspace=0.08, top=0.93, bottom=0.05)
    for column, (chrom, start, end) in enumerate(displayed):
        midpoint = (start + end) // 2
        zoom_start, zoom_end = midpoint - 1000, midpoint + 1000
        records = records_for_loci([(chrom, start, end)])[0]
        bands = read_cytobands(DATA / "mm10_cytoBandIdeo.txt.gz", chrom=chrom)
        parent_values_all = [
            query(DATA / filename, chrom, start, end) for _, filename, _ in SIGNALS
        ]
        child_values_all = [
            query(DATA / filename, chrom, zoom_start, zoom_end)
            for _, filename, _ in SIGNALS
        ]
        parent_max = max(1.0, max(float(values.max()) for values in parent_values_all) * 1.05)
        child_max = max(1.0, max(float(values.max()) for values in child_values_all) * 1.05)
        for row, (_, _filename, color) in enumerate(SIGNALS):
            parent_values = parent_values_all[row]
            child_values = child_values_all[row]
            parent_x = np.linspace(start, end, len(parent_values), endpoint=False)
            child_x = np.linspace(zoom_start, zoom_end, len(child_values), endpoint=False)
            parent_axis = axes[parent_row_start + row, column]
            child_axis = axes[child_row_start + row, column]
            plot_signal(parent_axis, parent_x, parent_values, color=color, alpha=1.0, linewidth=0)
            parent_axis.set_xlim(start, end)
            parent_axis.set_ylim(0, parent_max)
            parent_axis.axvspan(zoom_start, zoom_end, color=highlight_color, alpha=0.25, zorder=0)
            parent_axis.set_yticks([])
            plot_signal(child_axis, child_x, child_values, color=color, alpha=1.0, linewidth=0)
            child_axis.set_xlim(zoom_start, zoom_end)
            child_axis.set_ylim(0, child_max)
            child_axis.set_yticks([])
            if column == 0:
                parent_axis.set_ylabel(SIGNALS[row][0], fontsize=9)
                child_axis.set_ylabel(SIGNALS[row][0], fontsize=9)
            for axis in (parent_axis, child_axis):
                axis.tick_params(axis="x", bottom=False, labelbottom=False)
                style_axis(axis, hide_x=True)
                axis.spines[["top", "right", "left"]].set_color("#d0d0d0")
        if trap:
            trap_axis = axes[5, column]
            top_left = (zoom_start - start) / (end - start)
            top_right = (zoom_end - start) / (end - start)
            if zoom_side == "bottom":
                polygon_x = [top_left, top_right, 1.0, 0.0]
            else:
                polygon_x = [0.0, 1.0, top_right, top_left]
            trap_axis.fill(
                polygon_x, [0.9, 0.9, 0.05, 0.05],
                color="#888", alpha=0.55,
            )
            trap_axis.set_xlim(0, 1)
            trap_axis.set_ylim(0, 1)
            trap_axis.axis("off")
        draw_gene(axes[parent_gene_row, column], records, start, end)
        draw_gene(
            axes[child_gene_row, column], records, zoom_start, zoom_end,
            label_position="top" if zoom_side == "top" else "bottom",
        )
        axes[parent_gene_row, column].xaxis.set_visible(False)
        axes[child_gene_row, column].xaxis.set_visible(False)
        draw_cytoband(axes[11, column], bands, start, end)
        axes[11, column].xaxis.set_visible(False)
        if zoom_side == "bottom":
            axes[parent_row_start, column].set_title(
                f"{chrom}:{start:,}-{end:,}", fontsize=9, style="italic"
            )
        trap_axis = axes[5, column]
        trap_axis.text(
            0.5, 0.55, f"zoom: {chrom}:{zoom_start:,}-{zoom_end:,}",
            ha="center", va="center", fontsize=8, color="#777777",
        )
        axes[child_gene_row, column].set_xlabel(
            f"zoom {zoom_start:,}-{zoom_end:,}", fontsize=8
        )
    figure.suptitle(f"GenomeCanvas {name} parent/child zoom", fontsize=12, weight="bold")
    save_figure(figure, OUT, name, dpi=150)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panel_grid("p02_isoforms", LOCI14, gene_names=["Myh6", "Myh7", "Bcl2l2", "Pabpn1"])
    region_panel("p03_single_region", [(54_880_000, 54_930_000, "Region 1", "#66c2a5")])
    region_panel(
        "p03_multi_regions",
        [
            (54_842_000, 54_865_000, "Region A", "#66c2a5"),
            (54_875_000, 54_905_000, "Region B", "#fc8d62"),
            (54_920_000, 54_952_000, "Region C", "#8da0cb"),
            (54_965_000, 54_995_000, "Region D", "#e78ac3"),
        ],
    )
    region_panel(
        "p03_regions",
        [
            (54_842_000, 54_865_000, "Region A", "#66c2a5"),
            (54_875_000, 54_905_000, "Region B", "#fc8d62"),
            (54_920_000, 54_952_000, "Region C", "#8da0cb"),
            (54_965_000, 54_995_000, "Region D", "#e78ac3"),
        ],
    )
    panel_grid("p04_no_cyto", LOCI7, cyto=False)
    panel_grid("p05_with_cyto", LOCI7)
    panel_grid("p06_no_box", LOCI7, box=False)
    panel_grid("p07_with_box", LOCI7, box=True)
    panel_grid("p08_gene_bottom", LOCI7, cyto=False)
    panel_grid("p09_gene_top", LOCI7, cyto=False, gene_top=True)
    panel_grid("p10_ratio_08", LOCI7, aspect_ratio=0.8)
    panel_grid("p11_ratio_15", LOCI7, aspect_ratio=1.5)
    panel_grid("p12_yscale_gene", LOCI7, scale_mode="column")
    panel_grid("p13_yscale_track", LOCI7, scale_mode="row")
    panel_grid("p14_normal", LOCI7)
    panel_grid("p15_compact", LOCI7, compact=True)
    panel_grid("p16_default_colors", LOCI7)
    panel_grid("p17_custom_colors", LOCI7, colors=["#4c78a8", "#f58518", "#54a24b", "#e45756"])
    panel_grid("p18_flank_3k", [(chrom, start - 3000, end + 3000) for chrom, start, end in LOCI7])
    panel_grid(
        "p19_flank_10k", [(chrom, start - 10000, end + 10000) for chrom, start, end in LOCI7]
    )
    panel_grid("p20_label_bottom", LOCI7)
    panel_grid("p21_label_left", LOCI7)
    panel_grid("p22_highlights", P22_LOCI, box=True, highlights=P22_HIGHLIGHTS)
    trap_panel("p23_trap_default", gene_top=True)
    trap_panel("p24_trap_custom", gene_top=False, trap_color="#6baed6")
    zoom_panel("p25_zoom_f_bot", [P22_LOCI[0]], zoom_side="bottom")
    zoom_panel("p26_zoom_f_top", [P22_LOCI[1]], zoom_side="top")
    zoom_panel("p27_zoom_i_bot", [LOCI14[0]], zoom_side="bottom")
    zoom_panel("p28_zoom_i_top", [LOCI14[1]], trap=True, zoom_side="top")
    zoom_panel("p29_regions_zoom_bot", REGION_ZOOM_LOCI, zoom_side="bottom")
    zoom_panel("p30_regions_zoom_top", REGION_ZOOM_LOCI, trap=True, zoom_side="top")


if __name__ == "__main__":
    main()
