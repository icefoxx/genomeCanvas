"""Matplotlib heatmap and profile figures for enrichment matrices."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from ..render.plotting import save_figure
from .matrix import EnrichmentMatrix
from .plot_style import (
    HeatmapPlotStyle,
    ProfilePlotStyle,
    format_distance,
    mean_and_sem,
    scale_values,
)


def plot_enrichment_heatmap(
    matrix: EnrichmentMatrix,
    output: str | Path,
    *,
    zmin: float | None = None,
    zmax: float | None = None,
    show_profile: bool = True,
    height: float = 8.0,
    style: HeatmapPlotStyle | None = None,
    title: str | None = None,
    xlab: str | None = None,
    ylab: str | None = None,
    labels: list[str] | tuple[str, ...] | None = None,
    colorbar_label: str | None = None,
) -> Path:
    """Draw one column per sample: optional profile, heatmap, and colour key."""
    opts = style or HeatmapPlotStyle()
    if zmin is not None:
        opts.zmin = zmin
        opts.autoscale = False
    if zmax is not None:
        opts.zmax = zmax
        if opts.zmin is None:
            opts.zmin = 0.0
        opts.autoscale = False
    if not show_profile:
        opts.show_profile = False
    if height != 8.0:
        opts.height = height
    if title is not None:
        opts.title = title
    if xlab is not None:
        opts.xlab = xlab
    if ylab is not None:
        opts.ylab = ylab
    if labels is not None:
        opts.labels = labels
    if colorbar_label is not None:
        opts.colorbar_label = colorbar_label

    n_cols = max(1, len(matrix.sample_labels))
    sample_titles = _sample_titles(matrix, opts)
    positions = matrix.axis_positions()
    x0 = float(positions[0]) if positions.size else 0.0
    x1 = _right_edge(matrix, positions)
    n_rows = int(matrix.values.shape[0])
    width = max(3.6, opts.width_per_sample * n_cols + (0.7 if _show_groups(matrix, opts) else 0.0))
    ratios: list[float] = []
    if opts.show_profile:
        ratios.append(opts.profile_ratio)
    ratios.append(1.0)
    if opts.show_colorbar and opts.independent_scale:
        ratios.append(opts.colorbar_ratio)
    extra = 1 if _show_groups(matrix, opts) else 0
    figure = plt.figure(figsize=(width, opts.height), facecolor="white")
    grid = figure.add_gridspec(
        len(ratios),
        n_cols + extra,
        height_ratios=ratios,
        width_ratios=[1.0] * n_cols + ([0.12] if extra else []),
        hspace=0.22 if opts.show_colorbar and opts.independent_scale else 0.10,
        wspace=0.18,
        left=0.12 if extra else 0.10,
        right=0.88 if (opts.show_colorbar and not opts.independent_scale) else 0.97,
        top=0.88 if opts.title else 0.93,
        bottom=0.14,
    )
    profile_axes = []
    heat_axes = []
    bar_axes = []
    heat_row = 1 if opts.show_profile else 0
    images = []
    shared_norm: Normalize | None = None
    if opts.show_colorbar and not opts.independent_scale:
        shared_vmin, shared_vmax = _global_limits(matrix, opts)
        shared_norm = Normalize(vmin=shared_vmin, vmax=shared_vmax)

    plt.rcParams["font.family"] = opts.font_family
    for column, label in enumerate(matrix.sample_labels):
        block = scale_values(matrix.sample_block(column), opts.plot_scale)
        cmap = _colormap(
            opts.cmap or (matrix.colormaps[column] if column < len(matrix.colormaps) else "Blues")
        )
        vmin, vmax = _panel_limits(block, opts, shared_norm)
        if opts.show_profile:
            axis = figure.add_subplot(grid[0, column])
            profile_axes.append(axis)
            _draw_profile_panel(
                axis,
                block,
                positions,
                matrix,
                opts,
                column,
                x0,
                x1,
            )
            if column:
                axis.tick_params(labelleft=False)
            axis.set_title(sample_titles[column], fontsize=opts.sample_label_size, pad=6)
        heat = figure.add_subplot(
            grid[heat_row, column], sharex=profile_axes[column] if profile_axes else None
        )
        heat_axes.append(heat)
        image = heat.imshow(
            block,
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation=opts.interpolation,
            extent=(x0, x1, n_rows, 0.0),
            rasterized=True,
        )
        images.append(image)
        if opts.ln_h:
            for bound in matrix.group_boundaries[1:-1]:
                heat.axhline(bound, color="white", linewidth=1.0, zorder=3)
        if opts.ln_v:
            _vertical_guides(heat, matrix, x0, x1)
        heat.set_xlim(x0, x1)
        heat.tick_params(labelsize=opts.tick_size, length=3, color="#666666")
        if not opts.show_row_numbers:
            heat.set_yticks([])
        if column == 0:
            heat.set_ylabel(opts.ylab, fontsize=opts.axis_label_size)
        else:
            heat.set_ylabel("")
            heat.tick_params(labelleft=False)
        if not opts.show_profile:
            heat.set_title(sample_titles[column], fontsize=opts.sample_label_size, pad=6)
        _apply_distance_ticks(heat, x0, x1, opts.x_unit, matrix)
        heat.spines[["top", "right"]].set_visible(False)
        if opts.show_colorbar and opts.independent_scale:
            cax = figure.add_subplot(grid[heat_row + 1, column])
            bar_axes.append(cax)
            bar = figure.colorbar(image, cax=cax, orientation="horizontal")
            bar.ax.tick_params(labelsize=opts.legend_size, length=2)
            bar.set_label(
                opts.colorbar_label if column == n_cols // 2 else "", fontsize=opts.legend_size
            )
        heat.set_xlabel("")
    if heat_axes:
        xlabel = opts.xlab if opts.xlab is not None else _default_xlab(matrix)
        figure.text(0.5, 0.035, xlabel, ha="center", va="top", fontsize=opts.axis_label_size)
        for axis in heat_axes:
            for spine in axis.spines.values():
                spine.set_color("#555555")
    if extra:
        group_ax = figure.add_subplot(grid[heat_row, n_cols])
        _draw_group_strip(group_ax, matrix, opts)
    if opts.show_colorbar and not opts.independent_scale and images:
        bar = figure.colorbar(images[-1], ax=heat_axes, fraction=0.035, pad=0.02)
        bar.set_label(opts.colorbar_label, fontsize=opts.legend_size)
        bar.ax.tick_params(labelsize=opts.legend_size)
    if opts.title:
        figure.suptitle(opts.title, fontsize=opts.title_size, fontweight="bold", y=0.98)
    path = Path(output)
    save_figure(figure, path.parent, path.stem, dpi=opts.dpi, close=True)
    return path.with_suffix(".png")


def plot_profile_lines(
    matrix: EnrichmentMatrix,
    output: str | Path,
    *,
    title: str | None = None,
    style: ProfilePlotStyle | None = None,
    xlab: str | None = None,
    ylab: str | None = None,
    labels: list[str] | tuple[str, ...] | None = None,
) -> Path:
    """Draw one mean profile line per sample, optionally with SEM ribbons."""
    opts = style or ProfilePlotStyle()
    if title is not None:
        opts.title = title
    if xlab is not None:
        opts.xlab = xlab
    if ylab is not None:
        opts.ylab = ylab
    if labels is not None:
        opts.labels = labels
    plt.rcParams["font.family"] = opts.font_family
    figure, axis = plt.subplots(figsize=opts.figsize, facecolor="white")
    positions = matrix.axis_positions()
    palette = list(opts.colors) if opts.colors else _auto_colors(len(matrix.sample_labels))
    names = _sample_titles(matrix, opts)
    x0 = float(positions[0]) if positions.size else 0.0
    x1 = _right_edge(matrix, positions)
    for index, label in enumerate(matrix.sample_labels):
        block = scale_values(matrix.sample_block(index), opts.plot_scale)
        mean, sem = mean_and_sem(block)
        x = _bin_centers(matrix, positions)[: mean.size] if positions.size else np.arange(mean.size)
        color = matrix.profile_colors[index] if index < len(matrix.profile_colors) else "auto"
        if color == "auto":
            color = palette[index % len(palette)]
        if opts.error_estimates and np.any(sem > 0):
            axis.fill_between(
                x, mean - sem, mean + sem, color=color, alpha=opts.fill_alpha, linewidth=0
            )
        axis.plot(x, mean, linewidth=opts.linewidth, label=names[index], color=color)
    if opts.ln_v:
        _vertical_guides(axis, matrix, x0, x1)
    if opts.ln_h is not None:
        axis.axhline(opts.ln_h, color="#bbbbbb", linewidth=0.8, linestyle=":")
    axis.set_xlim(*(opts.xlim or (x0, x1)))
    if opts.ylim:
        axis.set_ylim(*opts.ylim)
    axis.set_xlabel(
        opts.xlab if opts.xlab is not None else _default_xlab(matrix), fontsize=opts.axis_label_size
    )
    axis.set_ylabel(opts.ylab, fontsize=opts.axis_label_size)
    heading = opts.title
    if heading is None:
        heading = f"n = {len(matrix.coords)} regions"
    axis.set_title(heading, fontsize=opts.title_size, pad=8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(labelsize=opts.tick_size, length=4)
    _apply_distance_ticks(axis, x0, x1, opts.x_unit, matrix)
    if opts.legend:
        axis.legend(
            loc=opts.legend_pos,
            frameon=False,
            fontsize=opts.legend_size,
            bbox_to_anchor=(1.02, 1.0) if opts.legend_pos in {"best", "upper left"} else None,
            borderaxespad=0.0,
        )
        figure.subplots_adjust(right=0.78)
    path = Path(output)
    save_figure(figure, path.parent, path.stem, dpi=opts.dpi, close=True)
    return path.with_suffix(".png")


def _right_edge(matrix: EnrichmentMatrix, positions: np.ndarray) -> float:  # type: ignore[type-arg]
    """End of the last bin, so the axis spans the full window (e.g. +3 kb, not +2.95 kb)."""
    if not positions.size:
        return 1.0
    if matrix.body:
        return float(matrix.body + matrix.downstream)
    step = float(positions[1] - positions[0]) if positions.size > 1 else float(matrix.bin_size)
    return float(positions[-1] + step)


def _bin_centers(matrix: EnrichmentMatrix, positions: np.ndarray) -> np.ndarray:  # type: ignore[type-arg]
    edges = np.append(positions, _right_edge(matrix, positions))
    return (edges[:-1] + edges[1:]) / 2


def _sample_titles(
    matrix: EnrichmentMatrix, opts: HeatmapPlotStyle | ProfilePlotStyle
) -> list[str]:
    defaults = [label.replace("\n", " ") for label in matrix.sample_labels]
    if not opts.labels:
        return defaults
    titles = list(opts.labels)
    if len(titles) < len(defaults):
        titles.extend(defaults[len(titles) :])
    return [str(item) for item in titles[: len(defaults)]]


def _show_groups(matrix: EnrichmentMatrix, opts: HeatmapPlotStyle) -> bool:
    return bool(opts.show_group_labels and len(matrix.group_labels) > 1)


def _default_xlab(matrix: EnrichmentMatrix) -> str:
    if matrix.body:
        return "relative position"
    point = {"center": "center", "TSS": "TSS", "TES": "TES"}.get(matrix.ref_point, "reference")
    return f"distance from {point} (bp)"


def _panel_limits(
    block: np.ndarray,  # type: ignore[type-arg]
    opts: HeatmapPlotStyle,
    shared: Normalize | None,
) -> tuple[float, float]:
    if shared is not None:
        return float(shared.vmin), float(shared.vmax)  # type: ignore[arg-type]
    finite = block[np.isfinite(block)]
    if opts.autoscale or (opts.zmin is None and opts.zmax is None):
        if finite.size == 0:
            return 0.0, 1.0
        low = float(np.nanpercentile(finite, 1))
        high = float(np.nanpercentile(finite, 99))
        if high <= low:
            high = low + 1e-6
        return max(low, 0.0) if opts.plot_scale == "linear" else low, high
    vmin = 0.0 if opts.zmin is None else float(opts.zmin)
    vmax = float(np.nanpercentile(finite, 99) or 1.0) if opts.zmax is None else float(opts.zmax)
    return vmin, max(vmax, vmin + 1e-9)


def _global_limits(matrix: EnrichmentMatrix, opts: HeatmapPlotStyle) -> tuple[float, float]:
    scaled = scale_values(matrix.values, opts.plot_scale)
    return _panel_limits(scaled, opts, None)


def _draw_profile_panel(
    axis,
    block: np.ndarray,  # type: ignore[type-arg]
    positions: np.ndarray,  # type: ignore[type-arg]
    matrix: EnrichmentMatrix,
    opts: HeatmapPlotStyle,
    column: int,
    x0: float,
    x1: float,
) -> None:
    color = matrix.profile_colors[column] if column < len(matrix.profile_colors) else "auto"
    palette = _auto_colors(max(1, len(matrix.group_labels)))
    for index, (left, right, group) in enumerate(
        zip(
            matrix.group_boundaries[:-1],
            matrix.group_boundaries[1:],
            matrix.group_labels,
            strict=True,
        )
    ):
        if right <= left:
            continue
        mean, sem = mean_and_sem(block[left:right])
        x = _bin_centers(matrix, positions)[: mean.size] if positions.size else np.arange(mean.size)
        line_color = None if color == "auto" else color
        if line_color is None:
            line_color = palette[index % len(palette)]
        if opts.error_estimates and np.any(sem > 0):
            axis.fill_between(x, mean - sem, mean + sem, color=line_color, alpha=0.18, linewidth=0)
        # A fixed sample colour cannot separate BED groups, so vary the dash pattern instead.
        dashes = "-" if color == "auto" else ("-", "--", ":", "-.")[index % 4]
        axis.plot(
            x, mean, linewidth=opts.linewidth, label=group, color=line_color, linestyle=dashes
        )
    axis.set_xlim(x0, x1)
    if opts.profile_ymin is not None or opts.profile_ymax is not None:
        ymin, ymax = axis.get_ylim()
        axis.set_ylim(
            opts.profile_ymin if opts.profile_ymin is not None else ymin,
            opts.profile_ymax if opts.profile_ymax is not None else ymax,
        )
    if opts.ln_v:
        _vertical_guides(axis, matrix, x0, x1)
    axis.spines[["top", "right", "bottom"]].set_visible(False)
    axis.tick_params(labelbottom=False, bottom=False, labelsize=opts.tick_size, length=3)
    if column == len(matrix.sample_labels) - 1 and len(matrix.group_labels) > 1:
        axis.legend(frameon=False, fontsize=opts.legend_size, loc="upper right")


def _vertical_guides(axis, matrix: EnrichmentMatrix, x0: float, x1: float) -> None:
    if x0 <= 0 <= x1:
        axis.axvline(0, color="#888888", linewidth=0.7, linestyle="--", zorder=2)
    if matrix.body and x0 <= matrix.body <= x1:
        axis.axvline(matrix.body, color="#888888", linewidth=0.7, linestyle="--", zorder=2)


def _apply_distance_ticks(axis, x0: float, x1: float, unit: str, matrix: EnrichmentMatrix) -> None:
    ticks = [x0, x1]
    if x0 < 0 < x1:
        ticks.append(0.0)
    if matrix.body and x0 < matrix.body < x1:
        ticks.append(float(matrix.body))
    ticks = sorted({round(item, 6) for item in ticks})
    axis.set_xticks(ticks)
    labels = []
    for item in ticks:
        if abs(item) < 1e-9:
            labels.append(
                {"center": "center", "TSS": "TSS", "TES": "TES"}.get(matrix.ref_point, "0")
            )
        elif matrix.body and abs(item - matrix.body) < 1e-6:
            labels.append("TES" if matrix.ref_point == "TSS" else "end")
        else:
            labels.append(format_distance(item, unit))  # type: ignore[arg-type]
    axis.set_xticklabels(labels)


def _draw_group_strip(axis, matrix: EnrichmentMatrix, opts: HeatmapPlotStyle) -> None:
    n_rows = max(matrix.values.shape[0], 1)
    palette = _auto_colors(len(matrix.group_labels))
    axis.set_xlim(0, 1)
    axis.set_ylim(n_rows, 0)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    for index, (left, right, group) in enumerate(
        zip(
            matrix.group_boundaries[:-1],
            matrix.group_boundaries[1:],
            matrix.group_labels,
            strict=True,
        )
    ):
        axis.axhspan(left, right, color=palette[index % len(palette)], alpha=0.35)
        mid = (left + right) / 2
        axis.text(
            0.5,
            mid,
            group,
            ha="center",
            va="center",
            fontsize=opts.group_label_size,
            rotation=90,
            clip_on=True,
        )


def _colormap(spec: str):
    if "," in spec:
        colors = [item.strip() for item in spec.split(",") if item.strip()]
        if colors:
            cmap = LinearSegmentedColormap.from_list("custom", colors)
            cmap.set_bad("#ffffff")
            return cmap
    try:
        cmap = plt.get_cmap(spec)
    except ValueError:
        cmap = plt.get_cmap("Blues")
    cmap = cmap.copy() if hasattr(cmap, "copy") else cmap
    if hasattr(cmap, "with_extremes"):
        return cmap.with_extremes(bad="#ffffff")
    if hasattr(cmap, "set_bad"):
        cmap.set_bad("#ffffff")
    return cmap


def _auto_colors(count: int) -> list[tuple[float, float, float]]:
    cmap = plt.get_cmap("tab10")
    return [cmap(index % 10)[:3] for index in range(count)]
