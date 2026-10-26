"""Publication figures for ProfileSeg (vector PDF/SVG + PNG)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import metrics as M  # noqa: E402
from .analysis import Analysis  # noqa: E402

RC = {
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
}


def _save(fig, out: Path, name: str, dpi: int) -> list[Path]:
    paths = []
    for ext in ("pdf", "svg", "png"):
        path = out / f"{name}.{ext}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def _colors(a: Analysis):
    pal = list(a.config.palette) + ["#0072B2", "#D55E00"]
    return pal[0], pal[1]


def _kb(x, _pos=None):
    return "anchor" if abs(x) < 1e-9 else f"{x / 1000:g} kb"


def make_figures(a: Analysis, out: Path) -> list[Path]:
    out.mkdir(parents=True, exist_ok=True)
    c_ref, c_alt = _colors(a)
    ref_name, alt_name = a.config.reference, a.config.alternative
    written: list[Path] = []
    with plt.rc_context({**RC, "font.family": [a.config.font, "DejaVu Sans"]}):
        rel = np.isin(np.arange(len(a.ref_idx) + len(a.alt_idx)), np.arange(len(a.ref_idx)))
        groups = ((ref_name, c_ref, rel), (alt_name, c_alt, ~rel))

        # A: per-replicate raw profiles + group means
        fig, ax = plt.subplots(figsize=(4.2, 3.0))
        x = M.centers(a.curves["raw_edges"])
        for name, col, mask in groups:
            reps = a.curves["raw_samples"][mask]
            for r in reps:
                ax.plot(x, r, color=col, lw=0.6, alpha=0.45)
            mean = reps.mean(axis=0)
            se = (
                reps.std(axis=0, ddof=1) / np.sqrt(len(reps))
                if len(reps) > 1
                else np.zeros_like(mean)
            )
            ax.fill_between(x, mean - se, mean + se, color=col, alpha=0.18, lw=0)
            ax.plot(x, mean, color=col, lw=1.8, label=f"{name} (n={len(reps)})")
        ax.axvline(0, color="#9aa5b1", lw=0.6, ls="--")
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.xaxis.set_major_formatter(_kb)
        ax.set(
            xlabel="distance from anchor",
            ylabel="mean signal (per region)",
            title="A · replicate profiles (thin) and condition means ± SE",
        )
        ax.legend(frameon=False)
        written += _save(fig, out, "A_replicate_profiles", a.config.dpi)

        # B: raw vs background-corrected, area-normalised
        fig, axes = plt.subplots(1, 2, figsize=(7.6, 2.9))
        for name, col, mask in groups:
            axes[0].plot(
                x, a.curves["raw_samples"][mask].mean(axis=0), color=col, lw=1.6, label=name
            )
        e = a.curves["edges"]
        xc, wc = M.centers(e), M.widths(e)
        for name, col, mask in groups:
            dens = a.curves["corrected_samples"][mask] / wc
            for d in dens:
                axes[1].plot(xc, d * 1000, color=col, lw=0.6, alpha=0.4)
            axes[1].plot(xc, dens.mean(axis=0) * 1000, color=col, lw=1.8, label=name)
        axes[0].set(title="B · raw", ylabel="mean signal", xlabel="distance from anchor")
        axes[1].set(
            title="background-corrected, area-normalised",
            ylabel="density (per kb)",
            xlabel="distance from anchor",
        )
        for axis in axes:
            axis.axvline(0, color="#9aa5b1", lw=0.6, ls="--")
            axis.xaxis.set_major_locator(plt.MaxNLocator(5))
            axis.xaxis.set_major_formatter(_kb)
        axes[1].axhline(0, color="#cbd2d9", lw=0.6)
        axes[0].legend(frameon=False)
        written += _save(fig, out, "B_raw_vs_normalised", a.config.dpi)

        # C: difference curve with bootstrap band (descriptive)
        fig, ax = plt.subplots(figsize=(4.2, 2.8))
        k = 1000 / wc
        ax.fill_between(
            xc,
            a.curves["diff_low"] * k,
            a.curves["diff_high"] * k,
            color="#52606d",
            alpha=0.2,
            lw=0,
            label="hierarchical bootstrap band",
        )
        ax.plot(xc, a.curves["diff"] * k, color="#1f2933", lw=1.6, label=f"{alt_name} − {ref_name}")
        ax.axhline(0, color="#9aa5b1", lw=0.8)
        ax.axvline(0, color="#9aa5b1", lw=0.6, ls="--")
        ax.xaxis.set_major_formatter(_kb)
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.set(
            xlabel="distance from anchor",
            ylabel="Δ density (per kb)",
            title="C · shape difference (descriptive; no local test in v1)",
        )
        ax.legend(frameon=False, loc="best")
        written += _save(fig, out, "C_difference_curve", a.config.dpi)

        # E: effect sizes, amplitude and shape in separate panels
        rows = {r["metric"]: r for r in a.effects}
        panels = [
            ("Amplitude (log2 FC)", [m for m in ("log2_auc", "log2_centre_auc") if m in rows]),
            (
                "Shape (bp)",
                [m for m in ("centroid", "sd_width", "iqr_width", "naive_sd_width") if m in rows],
            ),
            ("Asymmetry (Δ fraction)", [m for m in ("asymmetry",) if m in rows]),
        ]
        panels = [p for p in panels if p[1]]
        fig, axes = plt.subplots(
            1, len(panels), figsize=(3.6 * len(panels), 2.8), squeeze=False, layout="constrained"
        )
        for axis, (title, names) in zip(axes[0], panels, strict=True):
            for i, m in enumerate(names):
                r = rows[m]
                y = len(names) - 1 - i
                col = "#9aa5b1" if r["family"] == "naive" else "#1f2933"
                if np.isfinite(r["welch_low"]):
                    axis.plot(
                        [r["welch_low"], r["welch_high"]],
                        [y, y],
                        color=col,
                        lw=2.2,
                        solid_capstyle="round",
                    )
                axis.plot(
                    [r["boot_low"], r["boot_high"]],
                    [y - 0.18, y - 0.18],
                    color=col,
                    lw=1,
                    ls=(0, (2, 1.5)),
                )
                axis.plot(r["estimate"], y, "o", color=col, ms=5)
            axis.axvline(0, color="#d55e00", lw=0.8)
            axis.set_yticks(range(len(names)), [n.replace("_", " ") for n in reversed(names)])
            axis.set_title(title)
        floor = next(iter(rows.values()))["perm_floor"]
        fig.suptitle(
            f"E · {alt_name} vs {ref_name}: Welch CI (solid), bootstrap CI (dashed); "
            f"permutation p floor = {floor:.2g}",
            fontsize=9,
            y=1.04,
        )
        written += _save(fig, out, "E_effect_sizes", a.config.dpi)
    return written


def benchmark_figure(summary: list[dict], out: Path, dpi: int = 300) -> list[Path]:
    out.mkdir(parents=True, exist_ok=True)
    scenarios = list(dict.fromkeys(r["scenario"] for r in summary))
    metrics = list(dict.fromkeys(r["metric"] for r in summary))
    grid = np.full((len(scenarios), len(metrics)), np.nan)
    expected = np.zeros_like(grid)
    for r in summary:
        i, j = scenarios.index(r["scenario"]), metrics.index(r["metric"])
        grid[i, j] = r["welch_flag_rate"]
        expected[i, j] = r["expected"] != 0
    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(1.0 + 0.75 * len(metrics), 0.8 + 0.42 * len(scenarios)))
        im = ax.imshow(grid, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        for i in range(len(scenarios)):
            for j in range(len(metrics)):
                ax.text(
                    j,
                    i,
                    f"{grid[i, j]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white" if grid[i, j] > 0.6 else "#1f2933",
                )
                if expected[i, j]:
                    ax.add_patch(
                        plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#d55e00", lw=1.4)
                    )
        ax.set_xticks(
            range(len(metrics)), [m.replace("_", " ") for m in metrics], rotation=40, ha="right"
        )
        ax.set_yticks(range(len(scenarios)), [s.replace("_", " ") for s in scenarios])
        ax.set_title(
            "F · call rate (Welch CI excludes 0); orange = truly changes (noise-free effect)\n"
            "naive = raw shape incl. background, shown to expose the confound",
            fontsize=8.5,
        )
        fig.colorbar(im, ax=ax, fraction=0.03, label="fraction of runs")
        paths = _save(fig, out, "F_benchmark_call_rates", dpi)
    return paths


def regional_figure(a, regional, out: Path) -> list[Path]:
    """G: per-region volcano plots (moderated t) for the main amplitude and shape metrics."""
    out.mkdir(parents=True, exist_ok=True)
    names = [n for n in ("log2_signal_relative", "log2_enrichment", "shift") if n in regional.tests]
    with plt.rc_context({**RC, "font.family": [a.config.font, "DejaVu Sans"]}):
        fig, axes = plt.subplots(1, len(names), figsize=(3.0 * len(names), 2.8), squeeze=False, layout="constrained")
        for ax, name in zip(axes[0], names, strict=True):
            t = regional.tests[name]
            ok = np.isfinite(t["p"])
            colors = np.where(t["call"] == "Higher", "#d55e00", np.where(t["call"] == "Lower", "#0072b2", "#b8c0c8"))
            ax.scatter(t["estimate"][ok], -np.log10(np.maximum(t["p"][ok], 1e-300)), s=4, c=colors[ok], lw=0, alpha=0.7)
            ax.axvline(0, color="#9aa5b1", lw=0.6)
            up, down = int(np.sum(t["call"] == "Higher")), int(np.sum(t["call"] == "Lower"))
            ax.set(title=f"{name.replace('_', ' ')}\n{up} higher · {down} lower (FDR < {a.config.region_fdr})",
                   xlabel=f"{a.config.alternative} − {a.config.reference}", ylabel="−log10 p (moderated)")
        fig.suptitle("G · per-region tests (empirical-Bayes moderated t, BH FDR)", fontsize=9)
        return _save(fig, out, "G_region_volcano", a.config.dpi)
