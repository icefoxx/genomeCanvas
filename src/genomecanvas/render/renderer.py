"""Deterministic first rendering backend."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter, LinearLocator

from ..core.intervals import GenomicInterval
from ..interact.hits import HitRegistry
from ..io.bed import read_bed
from ..io.chromhmm import read_chromhmm
from ..io.cytoband import read_cytobands
from ..io.genes import collapse_transcripts, overlapping_transcripts
from ..io.gtf import read_gtf
from ..io.signal import query_signal
from ..io.vcf import read_vcf
from ..model.figure_spec import FigureSpec, TrackSpec
from ..model.overlays import HitTarget
from ..model.scale import ScaleSpec
from ..model.style import FigureStyle, font_from_style, track_height
from ..pipeline.plot_style import format_genomic_position
from ..model.track import SignalSlice
from ..provenance.manifest import build_manifest, write_manifest
from ..tracks.loop import LoopTrack
from .locus import (
    draw_chromhmm,
    draw_cytoband,
    draw_genes,
    draw_heatmap_rows,
    draw_highlights,
    draw_loops,
    draw_motifs,
    draw_overlay,
    draw_overlay_annotation,
    draw_peak_signal,
    draw_peaks,
    draw_signal,
    draw_variants,
    style_track_axis,
)

SIGNAL_TYPES = {"signal", "peak_signal", "derived_signal"}


def _limits(
    spec: FigureSpec, tracks: list[TrackSpec], values: dict[str, np.ndarray[Any, Any]]
) -> dict[str, tuple[float, float]]:
    """Resolve fixed, per-track, shared, robust, or symmetric y limits."""
    config = spec.scales
    mode = str(config.get("mode", "auto_highest"))
    result: dict[str, tuple[float, float]] = {}
    signal_tracks = [
        track for track in tracks
        if track.type in SIGNAL_TYPES or track.type in {"overlay", "difference"}
    ]
    if mode == "fixed":
        shared = ScaleSpec("fixed", config.get("ymin"), config.get("ymax")).resolve(
            np.array([0.0, 1.0])
        )
        return {track.id: shared for track in signal_tracks}
    groups: dict[str, list[float]] = {}
    for track in signal_tracks:
        array = values.get(track.id)
        if array is None:
            continue
        finite = array[np.isfinite(array)]
        group = str(track.scale_group or "__all_signal_tracks__")
        groups.setdefault(group, []).extend(float(item) for item in finite)
    for track in signal_tracks:
        style = track.style
        if "ymin" in style and "ymax" in style:
            result[track.id] = ScaleSpec(
                "fixed", style.get("ymin"), style.get("ymax")
            ).resolve(np.array([0.0, 1.0]))
        elif "ymin" in style or "ymax" in style:
            array = values.get(track.id, np.array([0.0, 1.0]))
            auto = ScaleSpec("per_track").resolve(array)
            ymin = float(style["ymin"]) if "ymin" in style else auto[0]
            ymax = float(style["ymax"]) if "ymax" in style else auto[1]
            if ymax <= ymin:
                ymax = ymin + 1.0
            result[track.id] = (ymin, ymax)
        elif mode in {"shared_group", "auto_highest", "shared_robust"}:
            group = str(track.scale_group or "__all_signal_tracks__")
            policy = (
                ScaleSpec("robust", quantile=float(config.get("quantile", 0.99)))
                if mode == "shared_robust"
                else ScaleSpec("per_track")
            )
            pool = groups.get(group, [0.0, 1.0])
            result[track.id] = policy.resolve(np.array(pool if pool else [0.0, 1.0]))
        else:
            result[track.id] = ScaleSpec(
                mode, quantile=float(config.get("quantile", 0.99))
            ).resolve(values.get(track.id, np.array([0.0, 1.0])))
    return result


def render(spec: FigureSpec, output: str | Path) -> Path:
    """Render aligned locus tracks to PDF/SVG/PNG and write hit targets."""
    spec = spec.resolved()
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not spec.region.is_resolved():
        raise ValueError("figure region must be resolved to chrom/start/end before rendering")
    chrom = str(spec.region.chrom)
    start = spec.region.start
    end = spec.region.end
    assert start is not None and end is not None
    view = GenomicInterval(chrom, start, end)
    style = spec.style or FigureStyle()
    plt.rcParams["font.family"] = style.font_family
    visible = [track for track in spec.tracks if track.visibility]
    overlay_members = {
        member for track in visible if track.type == "overlay" for member in track.members
    }
    visible = [track for track in visible if track.id not in overlay_members]
    height = spec.layout.height
    if height is None and spec.layout.aspect_ratio is not None:
        height = spec.layout.width / spec.layout.aspect_ratio
    ratios = [
        track_height(track.type, track.style.get("height", track.height))
        for track in visible
    ] or [1.0]
    height = height or 1.15 * max(1.0, sum(ratios))
    figure, axes = plt.subplots(
        max(1, len(visible)),
        1,
        sharex=False,
        squeeze=False,
        figsize=(spec.layout.width, height),
        facecolor=style.background,
        gridspec_kw={"height_ratios": ratios, "hspace": spec.layout.track_gap},
    )
    flat_axes = axes[:, 0]
    by_id = {track.id: track for track in spec.tracks}
    slices: dict[str, SignalSlice] = {}
    overlay_values: dict[str, np.ndarray[Any, Any]] = {}
    for track in spec.tracks:
        if track.type in SIGNAL_TYPES and track.source is not None:
            slices[track.id] = query_signal(
                track.source.path, chrom, start, end,
                pixels=min(2000, max(50, end - start)),
                input_type=track.input_type, engine=spec.bigwig_engine,
            )
        elif track.type in {"overlay", "difference", "heatmap"}:
            member_values = []
            for member_id in track.members:
                if member_id not in by_id:
                    raise ValueError(f"{track.type} {track.id!r} references unknown {member_id!r}")
                member = by_id[member_id]
                if member.source is None:
                    continue
                member_slice = query_signal(
                    member.source.path, chrom, start, end,
                    pixels=min(2000, max(50, end - start)),
                    input_type=member.input_type, engine=spec.bigwig_engine,
                )
                slices[member_id] = member_slice
                member_values.append(member_slice.values)
            if track.type == "difference" and len(member_values) == 2:
                left, right = member_values[0], member_values[1]
                operation = str(track.style.get("operation", "difference"))
                if operation == "ratio":
                    overlay_values[track.id] = left / np.maximum(right, 1e-9)
                elif operation == "log2_ratio":
                    ratio = np.maximum(left, 1e-9) / np.maximum(right, 1e-9)
                    overlay_values[track.id] = np.log2(ratio)
                else:
                    overlay_values[track.id] = left - right
                ref = slices[track.members[0]]
                slices[track.id] = SignalSlice(
                    chrom, start, end, ref.bin_starts, ref.bin_ends, overlay_values[track.id]
                )
            elif member_values and track.type == "overlay":
                overlay_values[track.id] = np.concatenate(member_values)
    display_values = {
        track_id: _display_values(
            slice_.values, spec.layout.y_value_transform, spec.layout.y_divisor
        )
        for track_id, slice_ in slices.items()
    }
    for track_id, values in overlay_values.items():
        display_values[track_id] = _display_values(
            values, spec.layout.y_value_transform, spec.layout.y_divisor
        )
    limits = _limits(spec, spec.tracks, display_values)
    hits = HitRegistry()
    if style.title:
        figure.suptitle(style.title, **style.title_font.matplotlib(style.font_family))
    def _label_too_tall(index, track):
        font = font_from_style(track.style, style.track_label_font, style.font_family)
        track_inches = height * 0.8 * ratios[index] / sum(ratios)
        return len(track.label or track.id) * font.size * 0.55 / 72 > track_inches

    # One orientation for all data tracks: mixing rotated and horizontal labels looks broken.
    horizontal_labels = any(
        _label_too_tall(i, t) for i, t in enumerate(visible) if ratios[i] >= 0.7
    )
    for axes_index, (axis, track) in enumerate(zip(flat_axes, visible, strict=False)):
        label_font = font_from_style(track.style, style.track_label_font, style.font_family)
        label_kwargs = label_font.matplotlib(style.font_family)
        if ratios[axes_index] < 0.7 or horizontal_labels:
            label_kwargs.update(rotation=0, ha="right", va="center")
        axis.set_ylabel(track.label or track.id, **label_kwargs)
        hits.add(
            HitTarget(
                id=f"track:{track.id}", kind="track", track_id=track.id,
                chrom=chrom, start=start, end=end, text=track.label or track.id,
                axes_index=axes_index, x=float(start), y=0.0,
            )
        )
        _draw_track(
            axis, track, spec, view, slices, display_values, limits, by_id, style, hits, axes_index,
        )
        if track.type != "cytoband":
            axis.set_xlim(start, end)
        if spec.highlights and track.type != "cytoband":
            draw_highlights(axis, spec.highlights, axis.get_ylim()[1])
        for overlay in spec.overlays:
            if overlay.track_id not in {None, track.id}:
                continue
            if overlay.chrom not in {None, chrom}:
                continue
            draw_overlay_annotation(
                axis, overlay, font=style.overlay_font, ylim=axis.get_ylim()
            )
            hits.add(
                HitTarget(
                    id=f"overlay:{overlay.id}", kind="overlay", overlay_id=overlay.id,
                    track_id=track.id, chrom=chrom,
                    start=int(overlay.start or start), end=int(overlay.end or overlay.start or end),
                    text=overlay.text, axes_index=axes_index,
                    x=float(overlay.x if overlay.x is not None else overlay.start or start),
                    y=float(overlay.y if overlay.y is not None else 0.0),
                )
            )
        style_track_axis(
            axis, show_frame=spec.layout.show_major_frame, tick=style.tick_font,
            spine_color=style.spine_color, grid=style.grid,
        )
        if not spec.layout.show_y_axis and track.type in SIGNAL_TYPES | {"overlay", "difference"}:
            axis.tick_params(axis="y", left=False, labelleft=False)
        if axes_index < len(flat_axes) - 1:
            axis.tick_params(axis="x", labelbottom=False, bottom=False)
    xlabel_axis = flat_axes[-1]
    if spec.layout.show_x_axis:
        span = end - start
        ticks = np.linspace(start, end, spec.layout.x_tick_count)
        xlabel_axis.set_xticks(ticks)
        xlabel_axis.set_xticklabels(
            [format_genomic_position(float(tick), span, spec.layout.x_tick_unit) for tick in ticks]
        )
        xlabel_axis.set_xlabel(
            spec.layout.x_label or f"{chrom}:{start:,}–{end:,}",
            **style.axis_label_font.matplotlib(style.font_family),
        )
    else:
        for axis in flat_axes:
            axis.tick_params(axis="x", labelbottom=False, bottom=False)
    figure.subplots_adjust(hspace=spec.layout.track_gap)
    if style.title:
        # Keep a fixed ~0.6 in title band instead of Matplotlib's 12% default.
        figure.subplots_adjust(top=1 - min(0.12, 0.6 / height))
    base = output_path.with_suffix("")
    figure.savefig(output_path, facecolor=style.background)
    figure.savefig(base.with_suffix(".pdf"), format="pdf", facecolor=style.background)
    figure.savefig(base.with_suffix(".svg"), format="svg", facecolor=style.background)
    plt.close(figure)
    spec.to_yaml(output_path.with_suffix(".spec.yaml"))
    write_manifest(build_manifest(spec), output_path.with_suffix(".manifest.json"))
    Path(str(base) + ".hits.json").write_text(
        json.dumps(hits.to_json(), indent=2) + "\n", encoding="utf-8"
    )
    return output_path


def _draw_track(
    axis: Any,
    track: TrackSpec,
    spec: FigureSpec,
    view: GenomicInterval,
    slices: dict[str, SignalSlice],
    display_values: dict[str, np.ndarray[Any, Any]],
    limits: dict[str, tuple[float, float]],
    by_id: dict[str, TrackSpec],
    style: FigureStyle,
    hits: HitRegistry,
    axes_index: int,
) -> None:
    chrom, start, end = view.chrom, view.start, view.end
    gene_font = font_from_style(track.style, style.gene_label_font, style.font_family)
    if track.type in SIGNAL_TYPES and track.id in slices:
        slice_ = slices[track.id]
        displayed = SignalSlice(
            slice_.chrom, slice_.start, slice_.end, slice_.bin_starts, slice_.bin_ends,
            display_values[track.id], slice_.units,
        )
        ymin, ymax = limits.get(track.id, (0.0, 1.0))
        if track.type == "peak_signal" and track.source is not None:
            draw_peak_signal(axis, read_bed(track.source.path, chrom=chrom, start=start, end=end))
        else:
            draw_signal(
                axis, displayed,
                color=str(track.style.get("color", "#3366aa")),
                ymin=ymin, ymax=ymax,
                fill=bool(track.style.get("fill", True)),
                fill_alpha=float(track.style.get("fill_alpha", style.signal_fill_alpha)),
                linewidth=float(track.style.get("linewidth", style.signal_linewidth)),
            )
        _configure_y(axis, spec, ymin, ymax)
    elif track.type == "difference" and track.id in slices:
        ymin, ymax = limits.get(track.id, (0.0, 1.0))
        draw_signal(
            axis, slices[track.id],
            color=str(track.style.get("color", "#8e44ad")),
            ymin=ymin, ymax=ymax, fill=bool(track.style.get("fill", True)),
        )
        _configure_y(axis, spec, ymin, ymax)
    elif track.type == "overlay":
        overlay_slices = []
        for member_id in track.members:
            member = by_id[member_id]
            if member_id not in slices:
                continue
            overlay_slices.append(
                (
                    member.label or member.id, slices[member_id],
                    str(member.style.get("color", track.style.get("color", "#3366aa"))),
                )
            )
        ymin, ymax = limits.get(track.id, (0.0, 1.0))
        draw_overlay(
            axis, overlay_slices, ymin=ymin, ymax=ymax,
            legend=style.legend_font,
            linewidth=float(track.style.get("linewidth", 1.2)),
        )
        _configure_y(axis, spec, ymin, ymax)
    elif track.type == "heatmap":
        rows = []
        for member_id in track.members:
            if member_id not in slices:
                continue
            member = by_id[member_id]
            color = str(member.style.get("color", "#3366aa"))
            rows.append((member.label or member.id, slices[member_id], color))
        draw_heatmap_rows(
            axis, rows, view=view, cmap=str(track.style.get("cmap", "Blues"))
        )
    elif track.type in {"peak", "annotation"} and track.source is not None:
        features = read_bed(track.source.path, chrom=chrom, start=start, end=end)
        draw_peaks(
            axis, features, color=str(track.style.get("color", "#2c3e50")),
            label_font=gene_font, show_labels=bool(track.style.get("show_labels", False)),
        )
        for feature in features:
            hits.add(
                HitTarget(
                    id=f"{track.type}:{track.id}:{feature.name}", kind=track.type,
                    track_id=track.id, chrom=chrom,
                    start=feature.interval.start, end=feature.interval.end,
                    text=feature.name, axes_index=axes_index,
                    x=float(feature.interval.start), y=0.5,
                )
            )
    elif track.type == "chromhmm" and track.source is not None:
        states = [
            (item.interval, item.state)
            for item in read_chromhmm(track.source.path)
            if item.interval.overlaps(view)
        ]
        raw_colors = track.style.get("colors")
        colors = (
            {str(key): str(value) for key, value in raw_colors.items()}
            if isinstance(raw_colors, dict) else None
        )
        draw_chromhmm(axis, states, colors=colors)
    elif track.type in {"gene", "isoform"}:
        source = track.source.path if track.source is not None else spec.gene_model
        if source is None:
            raise ValueError(
                f"{track.type} track {track.id!r} needs a GTF source or figure gene_model"
            )
        transcripts = overlapping_transcripts(read_gtf(source, chrom=chrom), chrom, start, end)
        collapse = spec.collapse_transcripts if track.type == "gene" else False
        if "collapse_transcripts" in track.style:
            collapse = bool(track.style.get("collapse_transcripts"))
        if collapse:
            transcripts = collapse_transcripts(transcripts)
        names = track.style.get("genes")
        if names:
            allowed = {str(name).casefold() for name in names}
            transcripts = [
                item for item in transcripts
                if item.gene_id.casefold() in allowed or (item.name or "").casefold() in allowed
            ]
        label_mode = str(track.style.get("label_mode", spec.layout.transcript_label))
        draw_genes(
            axis, transcripts, color=str(track.style.get("color", "#1b3a6b")),
            view=view, label_font=gene_font, label_mode=label_mode,
        )
        for transcript in transcripts:
            hits.add(
                HitTarget(
                    id=f"gene:{transcript.transcript_id}", kind="gene", track_id=track.id,
                    chrom=chrom, start=transcript.interval.start, end=transcript.interval.end,
                    text=transcript.name or transcript.transcript_id, axes_index=axes_index,
                    x=float(transcript.interval.start), y=0.0,
                    payload={"gene_id": transcript.gene_id},
                )
            )
    elif track.type == "cytoband":
        source = track.source.path if track.source is not None else None
        if source is None:
            raise ValueError(f"cytoband track {track.id!r} needs a source")
        bands = read_cytobands(source, chrom=chrom)
        draw_cytoband(axis, bands, view=view, marker_color=spec.layout.cytoband_marker_color)
    elif track.type == "loop" and track.source is not None:
        draw_loops(
            axis, LoopTrack(track.source.path).loops(), view=view,
            color=str(track.style.get("color", "#8e44ad")),
        )
    elif track.type == "variant" and track.source is not None:
        variants = read_vcf(track.source.path)
        draw_variants(
            axis, variants, view=view, color=str(track.style.get("color", "#c0392b")),
            label_font=gene_font,
        )
        for item in variants:
            pos = int(item["start"])
            hits.add(
                HitTarget(
                    id=f"variant:{item.get('id', pos)}", kind="variant", track_id=track.id,
                    chrom=str(item["chrom"]), start=pos, end=pos + 1, text=item.get("id"),
                    axes_index=axes_index, x=float(pos), y=0.7,
                )
            )
    elif track.type == "motif" and track.source is not None:
        features = read_bed(track.source.path, chrom=chrom, start=start, end=end)
        draw_motifs(
            axis, features, color=str(track.style.get("color", "#16a085")), label_font=gene_font
        )


def _configure_y(axis: Any, spec: FigureSpec, ymin: float, ymax: float) -> None:
    axis.set_ylim(ymin, ymax)
    ticks = np.linspace(ymin, ymax, spec.layout.y_tick_count)
    if ticks.size > 1:
        ticks = ticks[1:]
    axis.yaxis.set_major_locator(LinearLocator(spec.layout.y_tick_count))
    axis.set_yticks(list(ticks))
    axis.tick_params(
        axis="y",
        length=spec.layout.y_tick_length if spec.layout.show_y_ticks else 0,
        pad=spec.layout.y_tick_label_pad,
        labelsize=spec.style.tick_font.size,
        labelcolor=spec.style.tick_font.color,
    )
    axis.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value:.{spec.layout.y_tick_decimals}f}")
    )


def _display_values(
    values: np.ndarray[Any, Any], transform: str, divisor: float
) -> np.ndarray[Any, Any]:
    """Apply a renderer-only y-value display transform."""
    if transform == "divide":
        return values / divisor
    if transform == "log2":
        return np.asarray(np.log2(np.maximum(values, np.finfo(float).tiny)))
    if transform == "log10":
        return np.asarray(np.log10(np.maximum(values, np.finfo(float).tiny)))
    return values
