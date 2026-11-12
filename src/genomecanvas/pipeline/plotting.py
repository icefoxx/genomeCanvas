"""Matplotlib heatmap and profile figures for enrichment matrices."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
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
    label_pos = _label_pos(matrix, opts)
    extra = 1 if label_pos.startswith("strip") else 0
    strip_left = label_pos == "strip-left"
    shift = 1 if strip_left else 0  # grid column of the first sample
    figure = plt.figure(figsize=(width, opts.height), facecolor="white")
    grid = figure.add_gridspec(
        len(ratios),
        n_cols + extra,
        height_ratios=ratios,
        width_ratios=([0.12] if strip_left else [])
        + [1.0] * n_cols
        + ([0.12] if extra and not strip_left else []),
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
    blocks = [
        scale_values(matrix.sample_block(column), opts.plot_scale)
        for column in range(len(matrix.sample_labels))
    ]
    limits = heatmap_limits(blocks, opts, matrix.meta_groups)

    plt.rcParams["font.family"] = opts.font_family
    for column, label in enumerate(matrix.sample_labels):
        block = blocks[column]
        cmap = _colormap(
            opts.cmap or (matrix.colormaps[column] if column < len(matrix.colormaps) else "Blues")
        )
        vmin, vmax = limits[column]
        if opts.show_profile:
            axis = figure.add_subplot(grid[0, column + shift])
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
            axis.set_title(sample_titles[column], fontsize=opts.sample_label_size, pad=6)
        heat = figure.add_subplot(
            grid[heat_row, column + shift], sharex=profile_axes[column] if profile_axes else None
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
        if column == 0 and label_pos in {"text-left", "strip-left"}:
            heat.set_ylabel("")  # the group names take the place of the "regions" label
        elif column == 0:
            heat.set_ylabel(opts.ylab, fontsize=opts.axis_label_size)
        else:
            heat.set_ylabel("")
            heat.tick_params(labelleft=False)
        text_column = {"text-left": 0, "text-right": n_cols - 1}.get(label_pos)
        if column == text_column:
            # Name each region group (e.g. clusters C1…Ck) at its block, nudged apart if crowded.
            rows_pt = heat.get_position().height * opts.height * 72 / max(n_rows, 1)
            if label_pos == "text-right":
                heat = heat.twinx() if opts.show_row_numbers else heat
                heat.set_ylim(n_rows, 0)
                heat.yaxis.tick_right()
                heat.spines["right"].set_visible(False)
            heat.set_yticks(
                _group_label_rows(matrix, rows_pt, opts.axis_label_size),
                list(matrix.group_labels),
            )
            heat.tick_params(
                axis="y",
                length=0,
                labelsize=opts.axis_label_size,
                labelleft=label_pos == "text-left",
                labelright=label_pos == "text-right",
            )
            heat = heat_axes[-1]
        if not opts.show_profile:
            heat.set_title(sample_titles[column], fontsize=opts.sample_label_size, pad=6)
        _apply_distance_ticks(heat, x0, x1, opts.x_unit, matrix)
        heat.spines[["top", "right"]].set_visible(False)
        if opts.show_colorbar and opts.independent_scale:
            cax = figure.add_subplot(grid[heat_row + 1, column + shift])
            bar_axes.append(cax)
            bar = figure.colorbar(image, cax=cax, orientation="horizontal")
            bar.ax.tick_params(labelsize=opts.legend_size, length=2)
            bar.set_label(
                opts.colorbar_label if column == n_cols // 2 else "", fontsize=opts.legend_size
            )
        heat.set_xlabel("")
    if profile_axes:
        _apply_profile_ylims(profile_axes, opts, matrix.meta_groups)
    if heat_axes:
        xlabel = opts.xlab if opts.xlab is not None else _default_xlab(matrix)
        figure.text(0.5, 0.035, xlabel, ha="center", va="top", fontsize=opts.axis_label_size)
        for axis in heat_axes:
            for spine in axis.spines.values():
                spine.set_color("#555555")
    if extra:
        group_ax = figure.add_subplot(grid[heat_row, 0 if strip_left else n_cols])
        _draw_group_strip(group_ax, matrix, opts, "left" if strip_left else "right")
    if opts.show_colorbar and not opts.independent_scale and images:
        pad = 0.02 if label_pos != "text-right" else 0.08  # room for the group names
        bar = figure.colorbar(images[-1], ax=heat_axes, fraction=0.035, pad=pad)
        bar.set_label(opts.colorbar_label, fontsize=opts.legend_size)
        bar.ax.tick_params(labelsize=opts.legend_size)
    legend = None
    if profile_axes and len(matrix.group_labels) > 1 and opts.group_legend != "none":
        legend = _group_legend(figure, profile_axes, matrix, opts)
    if opts.title and legend is not None and opts.group_legend == "top":
        figure.canvas.draw()
        top = (legend.get_window_extent().y1 + 4 * figure.dpi / 72) / figure.bbox.height
        figure.suptitle(opts.title, fontsize=opts.title_size, fontweight="bold", y=top, va="bottom")
    elif opts.title:
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
    if opts.legend and opts.legend_pos == "top":
        legend = axis.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, 1.0),
            ncol=len(names) if len(names) <= 4 else -(-len(names) // 2),  # at most 2 rows
            frameon=False,
            fontsize=opts.legend_size,
            borderaxespad=0.4,
            handlelength=1.6,
            columnspacing=1.2,
        )
        _title_above(figure, axis, legend, heading, opts.title_size)
    elif opts.legend:
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
    return _label_pos(matrix, opts).startswith("strip")


def _default_xlab(matrix: EnrichmentMatrix) -> str:
    if matrix.body:
        return "relative position"
    point = {"center": "center", "TSS": "TSS", "TES": "TES"}.get(matrix.ref_point, "reference")
    return f"distance from {point} (bp)"


def _auto_limits(block: np.ndarray, opts: HeatmapPlotStyle) -> tuple[float, float]:  # type: ignore[type-arg]
    """1st-99th percentile colour range of one panel (floored at 0 on a linear scale)."""
    finite = block[np.isfinite(block)]
    if finite.size == 0:
        return 0.0, 1.0
    low = float(np.nanpercentile(finite, 1))
    high = float(np.nanpercentile(finite, 99))
    if opts.plot_scale == "linear":
        low = max(low, 0.0)
    if high <= low:
        high = low + 1e-6
    return low, high


def scale_groups(
    values: Sequence[float],
    mode: str = "auto",
    ratio: float = 2.0,
    labels: Sequence[str] | None = None,
) -> list[int]:
    """Assign panels to scale-sharing groups; panels in one group get one common range.

    ``auto`` sorts the panels by ``values`` (a range maximum or span) and starts a new group
    whenever a value exceeds ``ratio`` times the smallest value of the current group, so every
    group spans at most ``ratio``-fold. ``group`` shares within equal ``labels`` (metaGrp).
    Group ids are numbered in column order.
    """
    n = len(values)
    if mode == "shared":
        raw = [0] * n
    elif mode == "independent":
        raw = list(range(n))
    elif mode == "group":
        if not labels or len(labels) < n:
            raw = list(range(n))
        else:
            ids: dict[str, int] = {}
            raw = [ids.setdefault(str(labels[i]), len(ids)) for i in range(n)]
    elif mode == "auto":
        if ratio < 1:
            raise ValueError("scale ratio must be >= 1")
        raw = [0] * n
        finite = [i for i in range(n) if np.isfinite(values[i]) and values[i] > 0]
        ordered = sorted(finite, key=lambda k: values[k])
        # Divisive split: cut a run at its largest fold gap until each run spans <= ratio.
        pending, runs = ([ordered] if ordered else []), []
        while pending:
            run = pending.pop()
            if len(run) < 2 or values[run[-1]] <= values[run[0]] * ratio:
                runs.append(run)
                continue
            gaps = [values[run[k + 1]] / values[run[k]] for k in range(len(run) - 1)]
            cut = int(np.argmax(gaps)) + 1
            pending.extend([run[cut:], run[:cut]])
        for gid, run in enumerate(runs):
            for i in run:
                raw[i] = gid
        gid = len(runs) - 1
        for i in range(n):
            if i not in finite:
                gid += 1
                raw[i] = gid
    else:
        raise ValueError(f"unknown scale mode {mode!r}; use auto, shared, independent or group")
    order: dict[int, int] = {}
    return [order.setdefault(g, len(order)) for g in raw]


def _pool_limits(
    limits: list[tuple[float, float]], groups: list[int]
) -> list[tuple[float, float]]:
    pooled: dict[int, tuple[float, float]] = {}
    for (low, high), group in zip(limits, groups, strict=True):
        if group in pooled:
            plow, phigh = pooled[group]
            pooled[group] = (min(plow, low), max(phigh, high))
        else:
            pooled[group] = (low, high)
    return [pooled[group] for group in groups]


def heatmap_limits(
    blocks: list,
    opts: HeatmapPlotStyle,
    meta_groups: Sequence[str] | None = None,
) -> list[tuple[float, float]]:
    """Colour range for every sample column (zscale mode, then zmin/zmax overrides)."""
    natural = [_auto_limits(block, opts) for block in blocks]
    mode = opts.zscale if opts.independent_scale else "shared"  # one colour bar = one range
    groups = scale_groups([high for _, high in natural], mode, opts.zscale_ratio, meta_groups)
    pooled = _pool_limits(natural, groups)
    if opts.autoscale:
        return pooled
    final = []
    for index, (low, high) in enumerate(pooled):
        user_low = _per_sample(opts.zmin, index)
        user_high = _per_sample(opts.zmax, index)
        low = 0.0 if user_low is None else user_low
        high = high if user_high is None else user_high
        final.append((low, max(high, low + 1e-9)))
    return final


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
    if opts.ln_v:
        _vertical_guides(axis, matrix, x0, x1)
    axis.spines[["top", "right", "bottom"]].set_visible(False)
    axis.tick_params(labelbottom=False, bottom=False, labelsize=opts.tick_size, length=3)


def _title_above(figure, axis, legend, heading: str, size: float) -> None:
    """Re-pad the axis title so it sits just above a legend placed over the axes."""
    figure.canvas.draw()
    legend_top = legend.get_window_extent().y1
    axes_top = axis.get_window_extent().y1
    axis.set_title(heading, fontsize=size, pad=(legend_top - axes_top) * 72 / figure.dpi + 4)


def _label_pos(matrix: EnrichmentMatrix, opts: HeatmapPlotStyle) -> str:
    """Where group names go beside the heatmap; "none" for a single group."""
    if len(matrix.group_labels) < 2:
        return "none"
    if not opts.show_group_labels and opts.group_label_pos.startswith("strip"):
        return "text-left"
    return opts.group_label_pos


def _group_legend(figure, profile_axes, matrix: EnrichmentMatrix, opts: HeatmapPlotStyle):
    """One legend for the region groups / clusters, centred above the profile row."""
    from matplotlib.lines import Line2D

    fixed = any(c != "auto" for c in matrix.profile_colors)
    palette = _auto_colors(len(matrix.group_labels))
    names = list(opts.group_legend_labels or [])
    handles = [
        Line2D(
            [],
            [],
            linewidth=opts.linewidth,
            color="#555555" if fixed else palette[index % len(palette)],
            linestyle=("-", "--", ":", "-.")[index % 4] if fixed else "-",
            label=names[index] if index < len(names) and names[index] else group,
        )
        for index, group in enumerate(matrix.group_labels)
    ]
    figure.canvas.draw()
    if opts.group_legend == "right":
        last = profile_axes[-1].get_window_extent()
        return figure.legend(
            handles=handles,
            loc="center left",
            bbox_to_anchor=(
                (last.x1 + 8 * figure.dpi / 72) / figure.bbox.width,
                (last.y0 + last.y1) / 2 / figure.bbox.height,
            ),
            frameon=False,
            fontsize=opts.legend_size,
            handlelength=2.2,
            borderaxespad=0.0,
        )
    titles_top = max(a.title.get_window_extent().y1 for a in profile_axes)
    y = (titles_top + 2 * figure.dpi / 72) / figure.bbox.height
    return figure.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, y),
        ncol=min(len(handles), 8),
        frameon=False,
        fontsize=opts.legend_size,
        handlelength=2.2,
        columnspacing=1.4,
        borderaxespad=0.0,
    )


def _per_sample(value, index: int) -> float | None:
    """Resolve a scalar-or-per-sample y-limit setting for one panel."""
    if value is None:
        return None
    if isinstance(value, (int, float, np.floating, np.integer)):
        return float(value)
    items = list(value)
    if len(items) == 1:  # one value broadcasts to every panel
        index = 0
    if index >= len(items) or items[index] is None:
        return None
    return float(items[index])


def _apply_profile_ylims(
    axes: list, opts: HeatmapPlotStyle, meta_groups: Sequence[str] | None = None
) -> None:
    """Set profile-panel y limits: grouped by ``profile_yaxis``, then per-panel overrides.

    A panel shows tick labels when its limits differ from the panel to its left, so an
    unlabeled panel always reads off the nearest labeled axis.
    """
    natural = [tuple(float(v) for v in axis.get_ylim()) for axis in axes]
    groups = scale_groups(
        [high - low for low, high in natural], opts.profile_yaxis, opts.profile_share_ratio,
        meta_groups,
    )
    pooled = _pool_limits(natural, groups)  # type: ignore[arg-type]
    final: list[tuple[float, float]] = []
    for index, (axis, (low, high)) in enumerate(zip(axes, pooled, strict=True)):
        user_low = _per_sample(opts.profile_ymin, index)
        user_high = _per_sample(opts.profile_ymax, index)
        low = low if user_low is None else user_low
        high = high if user_high is None else user_high
        if high <= low:
            high = low + (abs(low) * 1e-3 or 1e-6)
        axis.set_ylim(low, high)
        final.append((low, high))
    for index, axis in enumerate(axes):
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        same_as_left = index > 0 and np.allclose(final[index], final[index - 1], rtol=1e-9, atol=0)
        axis.tick_params(labelleft=not same_as_left)


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


def _spread(centers: list[float], gap: float, low: float, high: float) -> list[float]:
    """Nudge sorted 1-D label positions apart to at least ``gap``, staying within [low, high]."""
    out = list(centers)
    for i in range(1, len(out)):
        out[i] = max(out[i], out[i - 1] + gap)
    overflow = out[-1] - high if out else 0.0
    if overflow > 0:
        out[-1] -= overflow
        for i in range(len(out) - 2, -1, -1):
            out[i] = min(out[i], out[i + 1] - gap)
    return [max(low, y) for y in out]


def _group_label_rows(matrix: EnrichmentMatrix, rows_pt: float, size: float) -> list[float]:
    """Row positions for one label per group: block centres, spread so labels do not collide."""
    bounds = matrix.group_boundaries
    centers = [(a + b) / 2 for a, b in zip(bounds[:-1], bounds[1:], strict=True)]
    gap = size * 1.25 / max(rows_pt, 1e-9)
    return _spread(centers, gap, 0.0, float(bounds[-1]))


def _draw_group_strip(
    axis, matrix: EnrichmentMatrix, opts: HeatmapPlotStyle, side: str = "right"
) -> None:
    n_rows = max(matrix.values.shape[0], 1)
    palette = _auto_colors(len(matrix.group_labels))
    axis.set_xlim(0, 1)
    axis.set_ylim(n_rows, 0)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    rows_pt = axis.get_position().height * axis.figure.get_figheight() * 72 / n_rows
    bounds = matrix.group_boundaries
    for index, (left, right) in enumerate(zip(bounds[:-1], bounds[1:], strict=True)):
        axis.axhspan(left, right, color=palette[index % len(palette)], alpha=0.35)
    fits = all(
        (b - a) * rows_pt >= opts.group_label_size * (0.6 * len(g) + 0.4)
        for a, b, g in zip(bounds[:-1], bounds[1:], matrix.group_labels, strict=True)
    )
    if fits:  # every label fits inside its block, written along the strip
        for (a, b), group in zip(
            zip(bounds[:-1], bounds[1:], strict=True), matrix.group_labels, strict=True
        ):
            axis.text(
                0.5,
                (a + b) / 2,
                group,
                ha="center",
                va="center",
                fontsize=opts.group_label_size,
                rotation=90,
            )
        return
    # Thin blocks: horizontal labels just outside the strip, nudged apart where crowded.
    x, ha = (1.15, "left") if side == "right" else (-0.15, "right")
    for row, group in zip(
        _group_label_rows(matrix, rows_pt, opts.group_label_size), matrix.group_labels, strict=True
    ):
        axis.text(
            x, row, group, ha=ha, va="center", fontsize=opts.group_label_size, clip_on=False
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
