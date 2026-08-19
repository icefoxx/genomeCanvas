"""Deterministic first rendering backend."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt  # type: ignore[import-untyped]
import numpy as np
from matplotlib.ticker import FuncFormatter, LinearLocator  # type: ignore[import-untyped]

from ..core.intervals import GenomicInterval
from ..io.bedgraph import read_bedgraph
from ..model.figure_spec import FigureSpec, TrackSpec
from ..model.scale import ScaleSpec
from ..provenance.manifest import build_manifest, write_manifest


def _limits(
    spec: FigureSpec, tracks: list[TrackSpec], values: dict[str, np.ndarray[Any, Any]]
) -> dict[str, tuple[float, float]]:
    """Resolve fixed, per-track, shared, robust, or symmetric y limits."""
    config = spec.scales
    mode = str(config.get("mode", "auto_highest"))
    result: dict[str, tuple[float, float]] = {}
    if mode == "fixed":
        shared = ScaleSpec("fixed", config.get("ymin"), config.get("ymax")).resolve(
            np.array([0.0, 1.0])
        )
        return {track.id: shared for track in tracks if track.type == "signal"}
    groups: dict[str, list[float]] = {}
    for track in tracks:
        if track.type != "signal":
            continue
        array = values[track.id]
        finite = array[np.isfinite(array)]
        group = str(track.scale_group or "__all_signal_tracks__")
        groups.setdefault(group, []).extend(float(item) for item in finite)
    for track in tracks:
        if track.type != "signal":
            continue
        style = track.style
        if "ymin" in style or "ymax" in style:
            result[track.id] = ScaleSpec(
                "fixed", style.get("ymin"), style.get("ymax")
            ).resolve(np.array([0.0, 1.0]))
        elif mode in {"shared_group", "auto_highest", "shared_robust"}:
            group = str(track.scale_group or "__all_signal_tracks__")
            policy = ScaleSpec("robust", quantile=float(config.get("quantile", 0.99))) \
                if mode == "shared_robust" else ScaleSpec("per_track")
            result[track.id] = policy.resolve(
                np.array(groups[group])
            )
        else:
            result[track.id] = ScaleSpec(
                mode, quantile=float(config.get("quantile", 0.99))
            ).resolve(values[track.id])
    return result


def render(spec: FigureSpec, output: str | Path) -> Path:
    """Render signal tracks and aligned basic tracks to PDF/SVG/PNG."""
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    visible = [track for track in spec.tracks if track.visibility]
    height = spec.layout.height
    if height is None and spec.layout.aspect_ratio is not None:
        height = spec.layout.width / spec.layout.aspect_ratio
    height = height or 1.5 * max(1, len(visible))
    figure, axes = plt.subplots(
        max(1, len(visible)), 1, sharex=True, squeeze=False,
        figsize=(spec.layout.width, height),
    )
    figure.subplots_adjust(hspace=spec.layout.track_gap)
    flat_axes = axes[:, 0]
    signal_values: dict[str, np.ndarray[Any, Any]] = {}
    records_by_track: dict[str, list[tuple[GenomicInterval, float]]] = {}
    for track in visible:
        if track.type == "signal" and track.source is not None:
            records = read_bedgraph(
                track.source.path, chrom=spec.region.chrom,
                start=spec.region.start, end=spec.region.end,
            )
            records_by_track[track.id] = records
            signal_values[track.id] = _display_values(
                np.array([value for _, value in records]), spec.layout.y_value_transform,
                spec.layout.y_divisor,
            )
    limits = _limits(spec, visible, signal_values)
    for axis, track in zip(flat_axes, visible, strict=False):
        axis.set_ylabel(track.label or track.id)
        if track.type == "signal" and track.source is not None:
            records = records_by_track[track.id]
            for interval, value in records:
                display_value = _display_values(
                    np.array([value]), spec.layout.y_value_transform, spec.layout.y_divisor
                )[0]
                axis.plot(
                    [interval.start, interval.end], [display_value, display_value],
                    color=track.style.get("color", "#3366aa"),
                )
            if track.id in limits:
                axis.set_ylim(*limits[track.id])
            axis.yaxis.set_major_locator(LinearLocator(spec.layout.y_tick_count))
            axis.tick_params(
                axis="y", length=spec.layout.y_tick_length if spec.layout.show_y_ticks else 0,
                pad=spec.layout.y_tick_label_pad,
            )
            axis.yaxis.set_major_formatter(
                FuncFormatter(
                    lambda value, _: f"{value:.{spec.layout.y_tick_decimals}f}"
                )
            )
            if not spec.layout.show_y_axis:
                axis.tick_params(
                    axis="y", left=False, right=False, labelleft=False, labelright=False
                )
                axis.yaxis.set_visible(False)
        axis.set_xlim(spec.region.start, spec.region.end)
        if not spec.layout.show_major_frame:
            for spine in axis.spines.values():
                spine.set_visible(False)
        axis.spines[["top", "right", "left"]].set_visible(False)
    flat_axes[-1].set_xlabel(spec.layout.x_label or f"{spec.region.chrom} (0-based half-open)")
    if not spec.layout.show_x_axis:
        for axis in flat_axes:
            axis.tick_params(axis="x", labelbottom=False, bottom=False)
    figure.tight_layout()
    base = output_path.with_suffix("")
    figure.savefig(output_path)
    figure.savefig(base.with_suffix(".pdf"), format="pdf")
    figure.savefig(base.with_suffix(".svg"), format="svg")
    plt.close(figure)
    spec.to_yaml(output_path.with_suffix(".spec.yaml"))
    write_manifest(build_manifest(spec), output_path.with_suffix(".manifest.json"))
    return output_path


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
