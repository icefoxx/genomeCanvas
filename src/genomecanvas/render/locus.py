"""Matplotlib drawing helpers for locus-style genomic tracks."""

from __future__ import annotations

from typing import Any

import numpy as np
from matplotlib.axes import Axes
from matplotlib.patches import FancyArrowPatch, Rectangle

from ..core.intervals import GenomicInterval
from ..io.chromhmm import chromhmm_color, ucsc_score_color
from ..io.cytoband import Cytoband, stain_color
from ..layout.interval_packing import pack_intervals
from ..model.feature import Feature
from ..model.overlays import OverlayAnnotation
from ..model.style import FontSpec
from ..model.track import SignalSlice
from ..model.transcript import Transcript


def style_track_axis(
    axis: Axes, *, show_frame: bool, tick: FontSpec, spine_color: str, grid: bool
) -> None:
    """Apply the shared locus-track spine policy."""
    axis.spines[["top", "right", "left"]].set_visible(False)
    for spine in axis.spines.values():
        spine.set_color(spine_color)
    if not show_frame:
        for spine in axis.spines.values():
            spine.set_visible(False)
    axis.tick_params(axis="y", length=3, labelsize=tick.size, labelcolor=tick.color)
    axis.tick_params(axis="x", labelsize=tick.size, labelcolor=tick.color)
    if grid:
        axis.grid(axis="y", color="#e8edf2", linewidth=0.5)


def draw_signal(
    axis: Axes,
    slice_: SignalSlice,
    *,
    color: str,
    ymin: float,
    ymax: float,
    fill: bool = True,
    fill_alpha: float = 0.85,
    linewidth: float = 0.6,
) -> None:
    """Draw binned signal as a filled step trace."""
    if slice_.values.size == 0:
        axis.set_ylim(ymin, ymax)
        return
    x = np.append(slice_.bin_starts, slice_.bin_ends[-1])
    y = np.append(slice_.values, slice_.values[-1])
    if fill:
        axis.fill_between(x, 0.0, y, step="post", color=color, alpha=fill_alpha, linewidth=0)
    axis.step(x, y, where="post", color=color, linewidth=linewidth)
    axis.set_ylim(ymin, ymax)


def draw_overlay(
    axis: Axes,
    slices: list[tuple[str, SignalSlice, str]],
    *,
    ymin: float,
    ymax: float,
    legend: FontSpec,
    linewidth: float = 1.2,
) -> None:
    """Draw multiple signals as lines on one axis."""
    for name, slice_, color in slices:
        if slice_.values.size == 0:
            continue
        x = (slice_.bin_starts + slice_.bin_ends) / 2
        axis.plot(x, slice_.values, color=color, linewidth=linewidth, label=name)
    axis.set_ylim(ymin, ymax)
    axis.legend(
        loc="upper right",
        frameon=False,
        fontsize=legend.size,
        prop={
            "family": legend.family or "DejaVu Sans",
            "style": legend.style,
            "weight": legend.weight,
        },
    )


def draw_peaks(
    axis: Axes,
    features: list[Feature],
    *,
    color: str = "#2c3e50",
    label_font: FontSpec | None = None,
    show_labels: bool = False,
) -> None:
    """Draw peak intervals as thick bars on a single row."""
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    axis.axhline(0.5, color="#e5e5e5", linewidth=4, zorder=0)
    for feature in features:
        axis.add_patch(
            Rectangle(
                (feature.interval.start, 0.15),
                max(1, feature.interval.length),
                0.7,
                facecolor=color,
                edgecolor="none",
            )
        )
        if show_labels and feature.name and label_font:
            axis.text(
                feature.interval.start, 0.92, feature.name,
                **label_font.matplotlib(), ha="left", va="top", clip_on=True,
            )


def draw_peak_signal(axis: Axes, features: list[Feature]) -> None:
    """Draw narrowPeak/broadPeak intervals with UCSC score shading."""
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    for feature in features:
        axis.add_patch(
            Rectangle(
                (feature.interval.start, 0.05),
                max(1, feature.interval.length),
                0.9,
                facecolor=ucsc_score_color(float(feature.score or 0.0)),
                edgecolor="none",
            )
        )


def draw_chromhmm(
    axis: Axes,
    states: list[tuple[GenomicInterval, str]],
    *,
    colors: dict[str, str] | None = None,
) -> None:
    """Draw chromHMM states as a single colored bar."""
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    for interval, state in states:
        axis.add_patch(
            Rectangle(
                (interval.start, 0.05),
                max(1, interval.length),
                0.9,
                facecolor=chromhmm_color(state, colors),
                edgecolor="none",
            )
        )


def draw_genes(
    axis: Axes,
    transcripts: list[Transcript],
    *,
    color: str = "#1b3a6b",
    view: GenomicInterval,
    label_font: FontSpec | None = None,
    label_mode: str = "name",
) -> None:
    """Draw packed transcript structures with strand chevrons."""
    if not transcripts:
        axis.set_ylim(0, 1)
        axis.set_yticks([])
        return
    rows = pack_intervals((item.interval.start, item.interval.end) for item in transcripts)
    row_of = {index: row for row, members in enumerate(rows) for index in members}
    n_rows = max(1, len(rows))
    axis.set_ylim(-0.2, n_rows + 0.2)
    axis.set_yticks([])
    font = label_font or FontSpec(size=6, style="italic", color="#555555")
    for index, transcript in enumerate(transcripts):
        row = row_of.get(index, 0)
        y = n_rows - row - 0.55
        start = max(transcript.interval.start, view.start)
        end = min(transcript.interval.end, view.end)
        if end <= start:
            continue
        axis.plot([start, end], [y, y], color=color, linewidth=0.8, solid_capstyle="butt")
        for exon in transcript.exons:
            left = max(exon.start, view.start)
            right = min(exon.end, view.end)
            if right <= left:
                continue
            axis.add_patch(
                Rectangle((left, y - 0.28), right - left, 0.56, facecolor=color, edgecolor="none")
            )
        strand = transcript.interval.strand
        if strand in {"+", "-"}:
            marker = ">" if strand == "+" else "<"
            axis.plot(start + (end - start) * 0.5, y, marker=marker, color=color, markersize=5)
        if label_mode == "none":
            continue
        if label_mode == "id":
            label = transcript.transcript_id
        else:
            label = transcript.name or transcript.gene_id or transcript.transcript_id
        axis.text(start, y + 0.38, label, **font.matplotlib(), clip_on=True)


def draw_cytoband(
    axis: Axes,
    bands: list[Cytoband],
    *,
    view: GenomicInterval,
    marker_color: str = "#e33d32",
) -> None:
    """Draw a chromosome ideogram with the current view highlighted."""
    if not bands:
        axis.set_ylim(0, 1)
        axis.set_yticks([])
        return
    chrom_end = max(band.end for band in bands)
    axis.set_xlim(0, chrom_end)
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    for band in bands:
        axis.add_patch(
            Rectangle(
                (band.start, 0.2),
                band.end - band.start,
                0.55,
                facecolor=stain_color(band.stain),
                edgecolor="#444444",
                linewidth=0.3,
            )
        )
    axis.add_patch(
        Rectangle(
            (view.start, 0.05),
            max(1, view.end - view.start),
            0.85,
            fill=False,
            edgecolor=marker_color,
            linewidth=1.6,
        )
    )
    axis.set_xticks([])


def draw_loops(
    axis: Axes,
    loops: list[dict[str, int | str]],
    *,
    view: GenomicInterval,
    color: str = "#8e44ad",
) -> None:
    """Draw BEDPE loops as arcs between anchor midpoints."""
    axis.set_ylim(0, 1.05)
    axis.set_yticks([])
    for item in loops:
        if str(item["chrom1"]) != view.chrom or str(item["chrom2"]) != view.chrom:
            continue
        left = (int(item["start1"]) + int(item["end1"])) / 2
        right = (int(item["start2"]) + int(item["end2"])) / 2
        if right < left:
            left, right = right, left
        if right < view.start or left > view.end:
            continue
        xs = np.linspace(left, right, 40)
        span = max(right - left, 1)
        height = min(0.95, 0.25 + 0.7 * (span / max(view.length, 1)))
        ys = height * (1 - ((xs - left) / span * 2 - 1) ** 2)
        axis.plot(xs, ys, color=color, linewidth=1.2)


def draw_variants(
    axis: Axes,
    variants: list[dict[str, str]],
    *,
    view: GenomicInterval,
    color: str = "#c0392b",
    label_font: FontSpec | None = None,
) -> None:
    """Draw SNVs/indels as lollipops."""
    axis.set_ylim(0, 1.15)
    axis.set_yticks([])
    font = label_font or FontSpec(size=6)
    for item in variants:
        if item.get("chrom") != view.chrom:
            continue
        pos = int(item.get("start") or item.get("pos") or 0)
        if not (view.start <= pos < view.end):
            continue
        axis.plot([pos, pos], [0.05, 0.7], color=color, linewidth=0.8)
        axis.plot(pos, 0.7, marker="o", color=color, markersize=5)
        label = item.get("id") or f"{item.get('ref', '')}>{item.get('alt', '')}"
        axis.text(pos, 0.78, label, **font.matplotlib(), ha="center", va="bottom", clip_on=True)


def draw_motifs(
    axis: Axes,
    features: list[Feature],
    *,
    color: str = "#16a085",
    label_font: FontSpec | None = None,
) -> None:
    """Draw motif hits as triangles with optional names."""
    axis.set_ylim(0, 1.2)
    axis.set_yticks([])
    font = label_font or FontSpec(size=6)
    for feature in features:
        mid = (feature.interval.start + feature.interval.end) / 2
        axis.plot(mid, 0.35, marker="^", color=color, markersize=8)
        if feature.name:
            axis.text(
                mid, 0.55, feature.name, **font.matplotlib(),
                ha="center", va="bottom", clip_on=True,
            )


def draw_heatmap_rows(
    axis: Axes,
    rows: list[tuple[str, SignalSlice, str]],
    *,
    view: GenomicInterval,
    cmap: str = "Blues",
) -> None:
    """Draw member signals as stacked heatmap rows."""
    if not rows:
        axis.set_yticks([])
        return
    stacked = np.vstack([
        np.nan_to_num(item.values, nan=0.0) for _, item, _ in rows
    ])
    axis.imshow(
        stacked,
        aspect="auto",
        cmap=cmap,
        interpolation="nearest",
        extent=(view.start, view.end, len(rows), 0),
    )
    axis.set_yticks([index + 0.5 for index in range(len(rows))])
    axis.set_yticklabels([name for name, _, _ in rows], fontsize=7)


def draw_highlights(axis: Axes, highlights: list[dict[str, Any]], ymax: float) -> None:
    """Overlay translucent highlight windows on a track."""
    for item in highlights:
        start = float(item.get("start", 0))
        end = float(item.get("end", 0))
        if end <= start:
            continue
        axis.axvspan(
            start,
            end,
            color=str(item.get("color", "#ffc41a")),
            alpha=float(item.get("alpha", 0.28)),
            zorder=0,
            linewidth=0,
        )


def draw_overlay_annotation(
    axis: Axes,
    item: OverlayAnnotation,
    *,
    font: FontSpec,
    ylim: tuple[float, float],
) -> None:
    """Draw one overlay note onto an axes."""
    if item.kind in {"highlight", "region_note"} and item.start is not None:
        if item.end is None:
            return
        axis.axvspan(item.start, item.end, color=item.fill, alpha=item.alpha, zorder=1, linewidth=0)
    x = item.x if item.x is not None else item.start
    if x is None:
        return
    y = item.y if item.y is not None else ylim[1]
    if item.kind == "arrow" and item.end is not None:
        axis.add_patch(
            FancyArrowPatch(
                (x, y), (item.end, y),
                arrowstyle="->", mutation_scale=10, color=item.color, linewidth=1.0,
            )
        )
    if item.text:
        kwargs = font.matplotlib()
        if item.font_size:
            kwargs["fontsize"] = item.font_size
        if item.font_style:
            kwargs["fontstyle"] = item.font_style
        if item.font_weight:
            kwargs["fontweight"] = item.font_weight
        kwargs["color"] = item.color
        axis.text(x, y, item.text, ha=item.ha, va=item.va, zorder=5, **kwargs)
