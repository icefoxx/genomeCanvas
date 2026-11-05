"""Figures for :mod:`genomecanvas.analysis.sample_qc`.

Labels follow the sample-count tiers: every sample is named up to 30 samples per panel, above
that only flagged samples are; heatmaps show names up to 60 samples and colour strips beyond.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import colormaps
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.patches import Patch

from .sample_qc import STATUS_ORDER, QCResult

STATUS_COLORS = {
    "mark_mismatch": "#9e1b32",
    "low_quality": "#e07b39",
    "swap": "#6a3d9a",
    "mislabel": "#d1495b",
    "mislabel_batch": "#f2a65a",
    "undetermined": "#c3c9d0",
    "ok": "#8fb8c9",
}
FLAGGED = {"mark_mismatch", "low_quality", "swap", "mislabel", "mislabel_batch"}
LABEL_ALL = 30
NAMES_MAX = 60


def palette(names) -> dict:
    colors = [colormaps["tab20"](i) for i in range(20)] + [
        colormaps["tab20b"](i) for i in range(20)
    ]
    return {name: colors[i % len(colors)] for i, name in enumerate(sorted(set(names)))}


def _short(name: str, width: int = 28) -> str:
    return name if len(name) <= width else name[: width - 1] + "…"


def draw_overview(result: QCResult, figure) -> None:
    """Mark-level PCA, mark confusion and status counts per mark on an existing figure."""
    axes = figure.subplots(1, 3, gridspec_kw={"width_ratios": [1.25, 1, 1]})
    colors = palette(result.marks)
    status = [row["status"] for row in result.rows]
    ax = axes[0]
    xy = result.overview_pca
    for mark, color in colors.items():
        pick = [i for i, m in enumerate(result.marks) if m == mark]
        ax.scatter(xy[pick, 0], xy[pick, 1], s=16, color=color, label=mark, linewidths=0)
    bad = [i for i, s in enumerate(status) if s == "mark_mismatch"]
    if bad:
        ax.scatter(
            xy[bad, 0], xy[bad, 1], s=60, facecolors="none", edgecolors="black", linewidths=1
        )
        for i in bad[:25]:
            ax.annotate(
                _short(result.samples[i]),
                xy[i],
                fontsize=6,
                xytext=(4, 4),
                textcoords="offset points",
            )
    v = result.overview_variance
    ax.set_xlabel(f"PC1 ({v[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({v[1]:.1f}%)")
    ax.set_title(f"All samples by mark (n = {len(result.samples)})", fontsize=10)
    ax.legend(
        fontsize=6, frameon=False, ncol=1 if len(colors) <= 14 else 2, loc="best", markerscale=1.2
    )

    ax = axes[1]
    labels = sorted(result.confusion)
    predicted = labels + sorted({p for c in result.confusion.values() for p in c} - set(labels))
    counts = np.array([[result.confusion[m].get(p, 0) for p in predicted] for m in labels], float)
    share = counts / np.maximum(counts.sum(axis=1, keepdims=True), 1)
    ax.imshow(share, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(predicted)), predicted, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    for (p, q), value in np.ndenumerate(counts):
        if value:
            ax.text(
                q,
                p,
                int(value),
                ha="center",
                va="center",
                fontsize=6,
                color="white" if share[p, q] > 0.6 else "#1f2933",
            )
    ax.set_xlabel("predicted mark (kNN vote)")
    ax.set_ylabel("labelled mark")
    ax.set_title("Mark label check", fontsize=10)

    ax = axes[2]
    order = [info["mark"] for info in result.mark_summary]
    left = np.zeros(len(order))
    for state in STATUS_ORDER:
        widths = np.array(
            [sum(1 for r in result.rows if r["mark"] == m and r["status"] == state) for m in order]
        )
        if widths.any():
            ax.barh(
                range(len(order)),
                widths,
                left=left,
                color=STATUS_COLORS[state],
                label=state,
                height=0.75,
            )
        left += widths
    ax.set_yticks(range(len(order)), order, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("samples")
    ax.set_title("Status per mark", fontsize=10)
    ax.legend(fontsize=6, frameon=False, loc="lower right")


def pca_facets(result: QCResult, path: Path, dpi: int = 200) -> list[Path]:
    marks = list(result.per_mark)
    if not marks:
        return []
    cols = min(4, len(marks))
    rows = int(np.ceil(len(marks) / cols))
    figure = Figure(figsize=(4.2 * cols, 3.7 * rows + 0.6), layout="constrained")
    axes = np.atleast_1d(figure.subplots(rows, cols)).ravel()
    colors = palette(result.groups)
    for ax, mark in zip(axes, marks, strict=False):
        _pca_panel(ax, result, mark, colors, label_all=False)
    for ax in axes[len(marks) :]:
        ax.set_visible(False)
    figure.suptitle("Within-mark PCA, coloured by group (flagged samples outlined and named)")
    out = []
    for ext in ("png", "pdf"):
        target = path.with_suffix("." + ext)
        figure.savefig(target, dpi=dpi)
        out.append(target)
    return out


def _pca_panel(ax, result: QCResult, mark: str, colors: dict, *, label_all: bool) -> None:
    info = result.per_mark[mark]
    idx, xy, v = info["index"], info["pca"], info["variance"]
    status = [result.rows[i]["status"] for i in idx]
    for a, i in enumerate(idx):
        flagged = status[a] in FLAGGED
        ax.scatter(
            xy[a, 0],
            xy[a, 1],
            s=34 if flagged else 22,
            color=colors[result.groups[i]],
            edgecolors=STATUS_COLORS[status[a]] if flagged else "white",
            linewidths=1.4 if flagged else 0.5,
            zorder=3 if flagged else 2,
        )
        if flagged or (label_all and len(idx) <= LABEL_ALL):
            ax.annotate(
                _short(result.samples[i], 22),
                xy[a],
                fontsize=5.5,
                xytext=(3, 3),
                textcoords="offset points",
                color=STATUS_COLORS[status[a]] if flagged else "#3e4c59",
            )
    ax.set_title(f"{mark} · n = {len(idx)}", fontsize=9)
    ax.set_xlabel(f"PC1 ({v[0]:.0f}%)", fontsize=7)
    ax.set_ylabel(f"PC2 ({v[1]:.0f}%)", fontsize=7)
    ax.tick_params(labelsize=6)


def mark_pages(result: QCResult, path: Path) -> Path:
    """One page per mark: group-ordered correlation heatmap, group x group medians, PCA."""
    group_colors = palette(result.groups)
    batch_colors = palette([b for b in result.batches if b])
    with PdfPages(path) as pdf:
        for mark, info in result.per_mark.items():
            idx, r = info["index"], info["r"]
            order = sorted(
                range(len(idx)),
                key=lambda a: (
                    result.groups[idx[a]],
                    result.batches[idx[a]],
                    result.samples[idx[a]],
                ),
            )
            figure = Figure(figsize=(17, 6.6), layout="constrained")
            grid = figure.add_gridspec(1, 3, width_ratios=[1.25, 1, 1])
            heat = grid[0].subgridspec(
                2, 2, width_ratios=[0.035, 1], height_ratios=[0.035, 1], wspace=0.01, hspace=0.01
            )
            ax = figure.add_subplot(heat[1, 1])
            shown = r[np.ix_(order, order)]
            low = (
                float(np.nanpercentile(shown[~np.eye(len(order), dtype=bool)], 2))
                if len(order) > 1
                else 0
            )
            image = ax.imshow(
                shown, cmap="viridis", vmin=low, vmax=1, aspect="auto", interpolation="nearest"
            )
            names = [result.samples[idx[a]] for a in order]
            flagged = [result.rows[idx[a]]["status"] in FLAGGED for a in order]
            if len(order) <= NAMES_MAX:
                size = max(3.5, min(7, 260 / len(order)))
                ax.yaxis.tick_right()
                ax.set_yticks(range(len(order)), [_short(n, 34) for n in names], fontsize=size)
                for tick, bad in zip(ax.get_yticklabels(), flagged, strict=True):
                    if bad:
                        tick.set_color(STATUS_COLORS["mislabel"])
                        tick.set_fontweight("bold")
            else:
                ax.set_yticks([])
            ax.set_xticks([])
            strip = [group_colors[result.groups[idx[a]]] for a in order]
            for where, shape in ((heat[0, 1], (1, -1)), (heat[1, 0], (-1, 1))):
                s_ax = figure.add_subplot(where)
                s_ax.imshow(np.array(strip).reshape(*shape, 4), aspect="auto")
                s_ax.set_xticks([])
                s_ax.set_yticks([])
                if shape[0] == 1:
                    top = s_ax
            has_batch = any(result.batches[i] for i in idx)
            if has_batch:
                b_ax = top.inset_axes([0, 1.25, 1, 1])
                b_ax.imshow(
                    np.array(
                        [batch_colors.get(result.batches[idx[a]], (1, 1, 1, 1)) for a in order]
                    ).reshape(1, -1, 4),
                    aspect="auto",
                )
                b_ax.set_axis_off()
            figure.colorbar(image, ax=ax, location="bottom", shrink=0.5, pad=0.02, label="r")
            strips = "strips: batch (top), group" if has_batch else "strip: group"
            top.set_title(
                f"{mark}: ordered by group ({strips})",
                fontsize=9.5,
                loc="left",
                pad=22 if has_batch else 6,
            )

            ax = figure.add_subplot(grid[1])
            groups, medians, pairs = result.group_matrices[mark]
            ax.imshow(np.ma.masked_invalid(medians), cmap="viridis", vmin=low, vmax=1)
            ax.set_xticks(range(len(groups)), groups, rotation=60, ha="right", fontsize=6.5)
            ax.set_yticks(range(len(groups)), groups, fontsize=6.5)
            if len(groups) <= 20:
                for (p, q), value in np.ndenumerate(medians):
                    text = "–" if np.isnan(value) else f"{value:.2f}"
                    ax.text(
                        q,
                        p,
                        text,
                        ha="center",
                        va="center",
                        fontsize=5,
                        color="white"
                        if not np.isnan(value) and value < (low + 1) / 2
                        else "#1f2933",
                    )
            ax.set_title(
                "Group × group median r (diagonal: within group; – = one member)", fontsize=10
            )

            ax = figure.add_subplot(grid[2])
            _pca_panel(ax, result, mark, group_colors, label_all=True)
            present = sorted({result.groups[i] for i in idx})
            ax.legend(
                handles=[Patch(color=group_colors[g], label=g) for g in present], fontsize=5.5,
                frameon=False, ncol=2, loc="best",
            )  # fmt: skip
            pdf.savefig(figure)
    return path
