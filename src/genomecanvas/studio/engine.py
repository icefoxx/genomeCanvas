"""Figure engine behind GenomeCanvas Studio.

Everything here is pure computation: callers pass already-authorised file records
(absolute paths plus display metadata) and an output directory. No HTTP, no database.
"""

from __future__ import annotations

import colorsys
import dataclasses
import hashlib
import json
import os
import pickle
import re
from collections import Counter
from dataclasses import replace
from importlib import metadata
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb  # noqa: E402

from ..analysis import diffpeak, pca_scores, volcano_plot  # noqa: E402
from ..io.bed import read_bed  # noqa: E402
from ..io.gtf import read_gtf  # noqa: E402
from ..io.signal import query_signal  # noqa: E402
from ..model.figure_spec import FigureSpec  # noqa: E402
from ..pipeline.cluster import cluster_labels  # noqa: E402
from ..pipeline.compare import peak_center_comparisons  # noqa: E402
from ..pipeline.correlation import correlation_heatmap, pair_grid  # noqa: E402
from ..pipeline.matrix import (  # noqa: E402
    compute_region_matrix,
    load_matrix,
    save_matrix,
    sort_matrix,
)
from ..pipeline.metadata import SampleMeta  # noqa: E402
from ..pipeline.plot_style import HeatmapPlotStyle, ProfilePlotStyle  # noqa: E402
from ..pipeline.plotting import plot_enrichment_heatmap, plot_profile_lines  # noqa: E402
from ..pipeline.renorm import renormalize  # noqa: E402
from ..render.renderer import render  # noqa: E402
from ..tracks.difference import compare_signal  # noqa: E402
from ..transform.background import transform_values  # noqa: E402

CATALOG = [
    ("Genomic tracks", ["Stacked signal tracks", "Overlaid signals", "Gene / transcript models"]),
    (
        "Enrichment & heatmaps",
        [
            "Region enrichment heatmap",
            "Mean signal profile",
            "Profile + heatmap",
            "Clustered heatmap",
        ],
    ),
    (
        "Sample relationships",
        ["Correlation heatmap", "PCA scatter plot", "PairGrid scatter matrix"],
    ),
    ("Differential analysis", ["Differential peak volcano", "Peak-center comparison heatmap"]),
    (
        "Signal transformations",
        ["Difference", "Ratio", "Log₂ ratio", "Local background correction"],
    ),
]
# Minimum inputs per plot, mirrored in the UI so it can explain what is missing before a run.
REQUIREMENTS = [
    [dict(signal=1), dict(signal=2), dict(genes=1)],
    [dict(signal=1, regions=1)] * 4,
    [dict(signal=2, regions=1)] * 3,
    [dict(signal=4, regions=1), dict(signal=1, regions=2)],
    [dict(signal=2)] * 4,
]
# fmt: off
# Curated, print-safe categorical palettes. Colours are cycled when samples outnumber them.
PALETTES = {
    "studio": ("Studio", ["#2A7F78", "#C0714F", "#4F6D8F", "#D9A441", "#8C5C8E", "#6B9E5A", "#B85C6E", "#5E8FA8"]),
    "okabe": ("Okabe–Ito (colour-blind safe)", ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#5D6B7A"]),
    "npg": ("Nature (NPG)", ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148"]),
    "lancet": ("Lancet", ["#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6"]),
    "nejm": ("NEJM", ["#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97"]),
    "tol": ("Paul Tol muted", ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"]),
    "morandi": ("Morandi", ["#7D8B8F", "#B5838D", "#8E9AAF", "#A3B18A", "#C9A66B", "#9A8C98", "#6D8A96", "#B07D62"]),
    "nord": ("Nord", ["#5E81AC", "#BF616A", "#A3BE8C", "#D08770", "#B48EAD", "#88C0D0", "#EBCB8B", "#4C566A"]),
    "tableau": ("Tableau 10", ["#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F"]),
    "ink": ("Ink (greyscale)", ["#1F2933", "#52606D", "#7B8794", "#9AA5B1", "#3E4C59", "#616E7C"]),
    # Backward-compatible names from earlier request.json files.
    "teal": ("Teal / Terracotta", ["#388B83", "#BB8063"]),
    "blue": ("Blue / Violet", ["#5087BC", "#9375AD"]),
}
COLORMAPS = [
    "auto", "Blues", "Reds", "Greens", "Purples", "Oranges", "Greys", "YlGnBu", "YlOrRd", "GnBu",
    "PuBuGn", "BuPu", "viridis", "magma", "inferno", "plasma", "cividis", "rocket",
    "coolwarm", "RdBu_r", "RdYlBu_r", "Spectral_r", "PiYG", "BrBG",
]
DIVERGING = ["RdBu_r", "coolwarm", "bwr", "PuOr_r", "PiYG", "BrBG", "Spectral_r", "RdYlBu_r", "viridis", "magma"]
ROCKET = ["#03051A", "#4C1D4B", "#A11A5B", "#E83F3F", "#F69C73", "#FAEBDD"]
THEMES = {
    "clean": {"axes.spines.top": False, "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out"},
    "classic": {"axes.spines.top": True, "axes.spines.right": True},
    "grid": {"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
             "grid.color": "#e6e9ed", "grid.linewidth": 0.6, "axes.axisbelow": True},
    "minimal": {"axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
                "ytick.left": False, "axes.grid": True, "grid.axis": "y", "grid.color": "#eceff2"},
}
FONTS = [
    "DejaVu Sans", "Arial", "Helvetica", "Helvetica Neue", "Avenir Next", "Inter", "Roboto",
    "Source Sans Pro", "Liberation Sans", "Gill Sans", "Futura", "Times New Roman", "Georgia",
    "Palatino", "DejaVu Serif", "Menlo",
]
# fmt: on
SUFFIXES = {
    ".bw": "signal",
    ".bigwig": "signal",
    ".bedgraph": "signal",
    ".bdg": "signal",
    ".bg": "signal",
    ".bed": "regions",
    ".narrowpeak": "regions",
    ".broadpeak": "regions",
    ".gtf": "genes",
}
# rcParams a client may override. Anything that touches files, external programs
# (text.usetex runs LaTeX) or the backend is refused.
RC_ALLOWED_PREFIXES = (
    "axes.",
    "xtick.",
    "ytick.",
    "grid.",
    "lines.",
    "font.",
    "legend.",
    "patch.",
    "hatch.",
    "boxplot.",
    "errorbar.",
    "scatter.",
    "image.",
    "contour.",
    "figure.facecolor",
    "figure.edgecolor",
    "figure.titlesize",
    "figure.titleweight",
    "mathtext.fontset",
    "savefig.facecolor",
    "savefig.edgecolor",
    "text.color",
    "polaraxes.",
)
MAX_REGIONS = 20000  # per-bin matrices (heatmaps, profiles)
MAX_SUMMARY_ROWS = 3_000_000  # sample-by-feature tables (correlation, PCA, volcano)
MAX_GENOME_WINDOWS = 15_000_000  # before zero filtering


def relationship_mode(options):
    mode = str(options.get("rel_mode") or "regions")
    if mode not in {"regions", "bins", "windows"}:
        raise ValueError("Summarise over must be regions, bins or windows")
    return mode


# ---------------------------------------------------------------- file inspection


def _extent(chroms: Counter, spans: dict) -> dict | None:
    if not chroms:
        return None
    chrom = chroms.most_common(1)[0][0]
    start, end = spans[chrom]
    return dict(chrom=chrom, start=int(start), end=int(end))


def _looks_like_chromhmm(path: Path, features) -> bool:
    name = path.name.lower()
    if "chromhmm" in name or "hmm" in name or "segway" in name:
        return True
    names = [f.name for f in features[:200] if f.name]
    return bool(names) and all(re.fullmatch(r"(E?\d{1,2})(_\w+)?", n) for n in names)


def inspect_file(path: str | Path, role: str) -> dict:
    """Cheap summary used by the UI for validation, locus fitting and chromosome hints."""
    path = Path(path)
    info: dict = dict(size=path.stat().st_size)
    suffix = path.suffix.lower()
    try:
        if suffix in (".bw", ".bigwig"):
            import pyBigWig

            handle = pyBigWig.open(str(path))
            if handle is None or not handle.isBigWig():
                raise ValueError("not a valid BigWig file")
            sizes = handle.chroms()
            info.update(
                kind="BigWig",
                chroms=sorted(sizes, key=lambda c: -sizes[c])[:60],
                n_chroms=len(sizes),
            )
            if info["size"] < 64 * 1024**2:
                for chrom in sizes:
                    intervals = handle.intervals(chrom)
                    if intervals:
                        info["locus"] = dict(
                            chrom=chrom, start=intervals[0][0], end=intervals[-1][1]
                        )
                        break
            handle.close()
        elif role == "signal" or suffix in (".bedgraph", ".bdg", ".bg"):
            chroms, spans, rows = Counter(), {}, 0
            with path.open(encoding="utf-8") as handle:
                for line in handle:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    fields = line.split()
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    chroms[chrom] += 1
                    lo, hi = spans.get(chrom, (start, end))
                    spans[chrom] = (min(lo, start), max(hi, end))
                    rows += 1
            info.update(
                kind="bedGraph",
                rows=rows,
                chroms=[c for c, _ in chroms.most_common(60)],
                locus=_extent(chroms, spans),
            )
        elif role == "genes" or suffix == ".gtf":
            transcripts = read_gtf(path)
            chroms, spans = Counter(), {}
            for tx in transcripts:
                iv = tx.interval
                chroms[iv.chrom] += 1
                lo, hi = spans.get(iv.chrom, (iv.start, iv.end))
                spans[iv.chrom] = (min(lo, iv.start), max(hi, iv.end))
            info.update(
                kind="GTF",
                rows=len(transcripts),
                genes=len({tx.name for tx in transcripts}),
                chroms=[c for c, _ in chroms.most_common(60)],
                locus=_extent(chroms, spans),
            )
        else:
            features = read_bed(path)
            chroms, spans = Counter(), {}
            for feature in features:
                iv = feature.interval
                chroms[iv.chrom] += 1
                lo, hi = spans.get(iv.chrom, (iv.start, iv.end))
                spans[iv.chrom] = (min(lo, iv.start), max(hi, iv.end))
            kind = {".narrowpeak": "narrowPeak", ".broadpeak": "broadPeak"}.get(suffix, "BED")
            if kind == "BED" and _looks_like_chromhmm(path, features):
                kind = "chromHMM"
            info.update(
                kind=kind,
                rows=len(features),
                chroms=[c for c, _ in chroms.most_common(60)],
                locus=_extent(chroms, spans),
                stranded=any(f.interval.strand for f in features[:500]),
            )
    except Exception as exc:  # noqa: BLE001 - reported to the user as a file warning
        info["error"] = f"Could not read file: {exc}"
    return info


def gene_search(paths, query, limit=25):
    """Case-insensitive gene / transcript search across GTF files."""
    query = query.strip().lower()
    if len(query) < 2:
        return []
    merged: dict[tuple, dict] = {}
    for path in paths:
        for tx in read_gtf(path):
            name = (tx.name or tx.gene_id or "").lower()
            if query not in name and query not in tx.transcript_id.lower():
                continue
            iv = tx.interval
            key = (tx.name, iv.chrom)
            if key in merged:
                hit = merged[key]
                hit["start"], hit["end"] = min(hit["start"], iv.start), max(hit["end"], iv.end)
            else:
                merged[key] = dict(
                    name=tx.name, chrom=iv.chrom, start=iv.start, end=iv.end, strand=iv.strand
                )
    hits = sorted(
        merged.values(), key=lambda h: (not h["name"].lower().startswith(query), h["name"])
    )
    return hits[:limit]


# ---------------------------------------------------------------- styling helpers


def boolish(value, default=True):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).lower() in {"1", "true", "yes", "on", "t", "y"}


def optional_float(value):
    if value in (None, "", "auto"):
        return None
    number_ = float(value)
    if not np.isfinite(number_):
        raise ValueError("Expected a finite number")
    return number_


def number(options, key, default, low=None, high=None):
    value = options.get(key)
    if value in (None, ""):
        return default
    value = float(value)
    if (
        not np.isfinite(value)
        or (low is not None and value < low)
        or (high is not None and value > high)
    ):
        raise ValueError(f"{key.replace('_', ' ').capitalize()} must be between {low} and {high}")
    return value


def integer(value, name, low, high):
    number_ = int(float(value))
    if not low <= number_ <= high:
        raise ValueError(f"{name} must be between {low} and {high}")
    return number_


def valid_color(value, fallback):
    text = str(value).strip()
    if re.fullmatch(r"[0-9a-fA-F]{3}|[0-9a-fA-F]{6}", text):
        text = "#" + text
    try:
        return to_hex(to_rgb(text))
    except (ValueError, TypeError):
        return fallback


def shade(color, factor):
    """Scale lightness: factor < 1 darkens, > 1 lightens."""
    h, lightness, s = colorsys.rgb_to_hls(*to_rgb(color))
    return to_hex(colorsys.hls_to_rgb(h, max(0, min(1, lightness * factor)), s))


def tint(color, amount):
    """Mix a colour with white (amount 0 = unchanged, 1 = white)."""
    r, g, b = to_rgb(color)
    return to_hex((r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount))


def palette_colors(options):
    raw = str(options.get("custom_palette") or "").replace(";", ",")
    custom = [valid_color(c, None) for c in raw.split(",") if c.strip()]
    custom = [c for c in custom if c]
    if custom:
        return custom
    return PALETTES.get(str(options.get("palette") or "studio"), PALETTES["studio"])[1]


def cmap_spec(spec, color):
    """Resolve a heatmap colour-scale name: auto = white → sample colour → deep shade."""
    spec = str(spec or "auto").strip()
    if spec == "auto":
        return f"#ffffff,{shade(color, 1.05)},{shade(color, 0.45)}"
    if spec == "rocket":
        return ",".join(reversed(ROCKET))
    return spec


def resolve_cmap(spec, fallback):
    spec = str(spec or fallback).strip()
    if spec == "rocket":
        return LinearSegmentedColormap.from_list("rocket", ROCKET)
    if "," in spec:
        return LinearSegmentedColormap.from_list(
            "custom", [c.strip() for c in spec.split(",") if c.strip()]
        )
    try:
        return plt.get_cmap(spec)
    except ValueError as exc:
        raise ValueError(f"Unknown colour scale: {spec}") from exc


def coerce(default, value):
    if isinstance(default, bool):
        return boolish(value)
    if isinstance(default, int) and not isinstance(value, bool):
        return int(value)
    if isinstance(default, float):
        return float(value)
    if isinstance(default, tuple) or isinstance(value, list):
        return tuple(value) if isinstance(value, (list, tuple)) else value
    return value


def apply_overrides(target, values, label):
    """Set any documented dataclass field from the Advanced JSON panel."""
    if not values:
        return target
    if not isinstance(values, dict):
        raise ValueError(f"Advanced '{label}' must be a JSON object")
    names = {f.name for f in dataclasses.fields(target)}
    for key, value in values.items():
        if key not in names:
            raise ValueError(f"Unknown {label} option '{key}'. Valid: {', '.join(sorted(names))}")
        setattr(target, key, coerce(getattr(target, key), value) if value is not None else None)
    return target


def advanced(payload):
    raw = payload.get("advanced") or {}
    if isinstance(raw, str):
        raw = json.loads(raw) if raw.strip() else {}
    if not isinstance(raw, dict):
        raise ValueError("Advanced overrides must be a JSON object")
    allowed = {"rc", "heatmap", "profile", "layout", "style", "track_style"}
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(
            f"Unknown advanced section(s): {', '.join(sorted(unknown))}. "
            f"Use {', '.join(sorted(allowed))}"
        )
    return raw


def rc_params(options, extra):
    rc = dict(THEMES.get(str(options.get("theme") or "clean"), THEMES["clean"]))
    family = str(options.get("font_family") or "").strip()
    if family:
        rc["font.family"] = [family, "DejaVu Sans"]
    rc["svg.fonttype"] = "none"  # keep text editable in Illustrator/Inkscape
    rc["pdf.fonttype"] = 42
    for key, value in (extra or {}).items():
        if key not in matplotlib.rcParams:
            raise ValueError(f"Unknown Matplotlib rcParam '{key}'")
        if not key.startswith(RC_ALLOWED_PREFIXES):
            raise ValueError(f"rcParam '{key}' cannot be changed in the studio")
        rc[key] = value
    return rc


def figure_size(options):
    """Requested (width, height) in inches; either may be None (automatic)."""
    scale = 1 / 25.4 if options.get("size_unit") == "mm" else 1.0
    width, height = (
        optional_float(options.get("fig_width")),
        optional_float(options.get("fig_height")),
    )
    return (width * scale if width else None, height * scale if height else None)


def plot_styles(options, title, *, show_profile, extra):
    scale = options.get("plot_scale", "linear")
    if scale not in {"linear", "log2", "zscore"}:
        scale = "linear"
    allowed_legend = {
        "best",
        "upper right",
        "upper left",
        "lower left",
        "lower right",
        "right",
        "center left",
        "center right",
        "lower center",
        "upper center",
        "center",
    }
    legend_pos = str(options.get("legend_pos", "best"))
    if legend_pos not in allowed_legend:
        legend_pos = "best"
    font = str(options.get("font_family") or "DejaVu Sans")
    x_unit = (
        options.get("x_unit") if options.get("x_unit") in {"auto", "bp", "kb", "Mb"} else "auto"
    )
    manual_y = boolish(options.get("ylim_manual", options.get("yauto", False)), False)
    ylim = (
        (optional_float(options.get("ymin")), optional_float(options.get("ymax")))
        if manual_y
        else None
    )
    if ylim is not None and None in ylim:
        raise ValueError("Set both Y-axis limits, or switch off manual Y range")
    width, height = figure_size(options)
    dpi = int(number(options, "dpi", 300, 72, 1200))
    heatmap = HeatmapPlotStyle(
        title=title,
        xlab=str(options.get("xlab") or "") or None,
        ylab=str(options.get("ylab") or "regions"),
        colorbar_label=str(options.get("colorbar_label") or "signal"),
        show_profile=show_profile,
        show_colorbar=boolish(options.get("show_colorbar", True)),
        independent_scale=boolish(options.get("indi", True)),
        autoscale=boolish(options.get("autoscale", True)),
        zmin=optional_float(options.get("zmin")),
        zmax=optional_float(options.get("zmax")),
        plot_scale=scale,
        ln_v=boolish(options.get("ln_v", True)),
        error_estimates=boolish(options.get("error_estimates", True)),
        interpolation=str(options.get("interpolation") or "bilinear"),
        font_family=font,
        title_size=number(options, "title_size", 12, 4, 60),
        axis_label_size=number(options, "axis_label_size", 10, 4, 48),
        tick_size=number(options, "tick_size", 8, 3, 40),
        sample_label_size=number(options, "sample_label_size", 10, 3, 48),
        legend_size=number(options, "legend_size", 8, 3, 40),
        height=height or number(options, "height", 6, 2, 40),
        width_per_sample=number(options, "width_per_sample", 2.7, 0.8, 12),
        linewidth=number(options, "linewidth", 1.5, 0.1, 10),
        show_group_labels=boolish(options.get("show_group_labels", True)),
        profile_ymin=ylim[0] if ylim else None,
        profile_ymax=ylim[1] if ylim else None,
        x_unit=x_unit,
        cmap=None,
        dpi=dpi,
    )
    profile = ProfilePlotStyle(
        title=title,
        xlab=str(options.get("xlab") or "") or None,
        ylab=str(options.get("ylab") or "mean signal"),
        legend=boolish(options.get("legend", True)),
        legend_pos=legend_pos,  # type: ignore[arg-type]
        error_estimates=boolish(options.get("error_estimates", True)),
        ln_v=boolish(options.get("ln_v", True)),
        ln_h=optional_float(options.get("ln_h_pos"))
        if boolish(options.get("ln_h", False), False)
        else None,
        plot_scale=scale,
        ylim=ylim,
        font_family=font,
        title_size=number(options, "title_size", 12, 4, 60),
        axis_label_size=number(options, "axis_label_size", 11, 4, 48),
        tick_size=number(options, "tick_size", 9, 3, 40),
        legend_size=number(options, "legend_size", 9, 3, 40),
        linewidth=number(options, "linewidth", 1.8, 0.1, 10),
        fill_alpha=number(options, "fill_alpha", 0.18, 0, 1),
        x_unit=x_unit,
        figsize=(width or 7.2, height or 4.2),
        dpi=dpi,
    )
    apply_overrides(heatmap, extra.get("heatmap"), "heatmap")
    apply_overrides(profile, extra.get("profile"), "profile")
    return heatmap, profile


def parse_highlights(text, default_color):
    """Lines like `31136000-31138000 #f5c542 0.25` or `chr6:31,136,000-31,138,000`."""
    items = []
    for raw in str(text or "").replace(";", "\n").splitlines():
        line = raw.strip()
        if not line:
            continue
        match = re.match(
            r"^(?:[\w.]+:)?([\d,_]+)\s*[-–]\s*([\d,_]+)(?:\s+(\S+))?(?:\s+([\d.]+))?$", line
        )
        if not match:
            raise ValueError(f"Highlight '{line}' should look like 'start-end [#color] [alpha]'")
        start, end = (int(re.sub(r"[,_]", "", g)) for g in match.group(1, 2))
        items.append(
            dict(
                start=start,
                end=end,
                color=valid_color(match.group(3) or default_color, default_color),
                alpha=float(match.group(4) or 0.22),
            )
        )
    return items


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def software_versions():
    versions = {}
    for name in ("genomecanvas", "matplotlib", "numpy", "pyBigWig", "scipy", "scikit-learn"):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    return versions


def fit_to_page(fig, pad=0.04, rounds=8):
    """Shrink all axes so labels and titles fit inside a fixed page size.

    Exact-size export cannot crop to content, so overflow measured from the tight bounding
    box is reclaimed by scaling axes positions. Text does not scale, hence a few rounds.
    """
    fig.set_layout_engine(None)
    width, height = fig.get_size_inches()
    for _ in range(rounds):
        fig.canvas.draw()
        box = fig.get_tightbbox(fig.canvas.get_renderer())
        left, bottom = max(0.0, -box.x0) + pad, max(0.0, -box.y0) + pad
        right, top = max(0.0, box.x1 - width) + pad, max(0.0, box.y1 - height) + pad
        if max(left, bottom, right, top) <= pad + 1e-3:
            return
        sx = max(0.3, (width - left - right) / width)
        sy = max(0.3, (height - bottom - top) / height)
        for axis in fig.axes:
            p = axis.get_position()
            axis.set_position(
                [left / width + p.x0 * sx, bottom / height + p.y0 * sy, p.width * sx, p.height * sy]
            )


# ---------------------------------------------------------------- cache


def _fingerprint(path):
    stat = os.stat(path)
    return [str(path), stat.st_size, stat.st_mtime_ns]


class Cache:
    """Content-keyed cache for region matrices and summaries (restyling is then instant)."""

    def __init__(self, directory: Path | None):
        self.directory = directory
        if directory:
            directory.mkdir(parents=True, exist_ok=True)

    def key(self, kind, **parts):
        blob = json.dumps(parts, sort_keys=True, default=str).encode()
        return f"{kind}-{hashlib.sha256(blob).hexdigest()[:32]}"

    def path(self, key, suffix):
        return self.directory / f"{key}{suffix}" if self.directory else None

    def touch(self, path):
        try:
            os.utime(path)
        except OSError:
            pass


# ---------------------------------------------------------------- plotting


def run_plot(payload, records, destination, cache_dir=None):
    """Render one figure. `records` are authorised inputs with absolute `path`s."""
    destination = Path(destination)
    category = integer(payload.get("category", 1), "Category", 0, 4)
    plot = integer(payload.get("plot", 0), "Plot", 0, len(CATALOG[category][1]) - 1)
    if not records:
        raise ValueError("Select at least one input file")
    options = payload.get("settings", {}) or {}
    if not isinstance(options, dict):
        raise ValueError("Settings must be a JSON object")
    options = {**options, **({"palette": payload["palette"]} if payload.get("palette") else {})}
    extra = advanced(payload)
    colors = palette_colors(options)
    records = [dict(r) for r in records]
    for record in records:
        if record["role"] not in ("signal", "regions", "genes", "blacklist"):
            raise ValueError(f"Unknown role for {record['name']}")
    signals = [r for r in records if r["role"] == "signal"]
    regions = [r for r in records if r["role"] == "regions"]
    genes = [r for r in records if r["role"] == "genes"]
    blacklist = [r for r in records if r["role"] == "blacklist"]
    for index, record in enumerate(signals):
        record["color"] = valid_color(record.get("color") or "", colors[index % len(colors)])
    for index, record in enumerate(regions):
        record["color"] = valid_color(
            record.get("color") or "", colors[(len(signals) + index) % len(colors)]
        )
    if len(signals) > 16:
        raise ValueError("At most 16 signals per figure")
    flank = integer(options.get("flank", 500) or 500, "Flank", 20, 1_000_000)
    binsize = integer(options.get("binsize", 20) or 20, "Bin size", 1, 100_000)
    if flank * 2 // binsize > 4000:
        raise ValueError("Choose a larger bin size (maximum 4000 bins per window)")
    start = integer(options.get("start", 0) or 0, "Start", 0, 10**10)
    end = integer(options.get("end", 100000) or 100000, "End", start + 1, start + 20_000_000)
    chrom = str(options.get("chrom") or "chr1").strip()
    title = str(payload.get("title") or options.get("title") or CATALOG[category][1][plot])[:200]
    dpi = int(number(options, "dpi", 300, 72, 1200))
    transparent = boolish(options.get("transparent", False), False)
    fig_size = figure_size(options)
    exact = boolish(options.get("exact_size", False), False) and all(fig_size)
    output = destination / "figure.png"
    sources = {f"{r['label']}_{i + 1}": r["path"] for i, r in enumerate(signals)}
    names = list(sources)
    need = dict(REQUIREMENTS[category][plot])
    if category == 2 and relationship_mode(options) != "regions":
        need.pop("regions", None)  # genome bins / sliding windows need only signals
    have = dict(signal=len(signals), regions=len(regions), genes=len(genes))
    if any(have[k] < n for k, n in need.items()):
        parts = [f"{n} {k} file{'s' if n > 1 else ''}" for k, n in need.items()]
        raise ValueError(f"{CATALOG[category][1][plot]} needs at least " + " and ".join(parts))
    warnings: list[str] = []
    context = dict(
        category=category,
        plot=plot,
        options=options,
        extra=extra,
        title=title,
        signals=signals,
        regions=regions,
        genes=genes,
        blacklist=blacklist,
        names=names,
        sources=sources,
        colors=colors,
        chrom=chrom,
        start=start,
        end=end,
        flank=flank,
        binsize=binsize,
        destination=destination,
        output=output,
        warnings=warnings,
        fig_size=fig_size,
        cache=Cache(Path(cache_dir) if cache_dir else None),
        notes={},
        linewidth=number(options, "linewidth", 1.5, 0.1, 10),
        fill_alpha=number(options, "fill_alpha", 0.25, 0, 1),
    )

    # Capture helper-generated figures before they close to add title and SVG uniformly.
    def export_current(fig=None):
        fig = fig if hasattr(fig, "savefig") else plt.gcf()
        if exact:
            fig.set_size_inches(*fig_size)
            fit_to_page(fig)
        elif category in (2, 3, 4) and any(fig_size):
            width, height = fig.get_size_inches()
            fig.set_size_inches(fig_size[0] or width, fig_size[1] or height)
        has_axis_title = any(axis.get_title() for axis in fig.axes)
        if category == 2 and plot == 1 and fig.axes:
            fig.axes[0].margins(x=0.18, y=0.14)  # room for point labels
        if fig._suptitle is None and title and not has_axis_title:
            size = number(options, "title_size", 12, 4, 60)
            plain = [a for a in fig.axes if a.get_label() != "<colorbar>"]
            if len(plain) == 1:
                plain[0].set_title(title, fontsize=size, pad=10)  # layout-aware, never overlaps
            else:
                fig.suptitle(title, fontsize=size)
        for ext in ("svg", "png", "pdf"):
            fig.savefig(
                destination / f"figure.{ext}",
                dpi=dpi,
                transparent=transparent,
                bbox_inches=None if exact else "tight",
                facecolor="none" if transparent else "white",
            )

    original_close = plt.close

    def close_and_export(fig=None):
        if plt.get_fignums():
            export_current(fig)
        original_close(fig)

    plt.close = close_and_export
    try:
        with plt.rc_context(rc_params(options, extra.get("rc"))):
            if category == 0:
                draw_locus(context)
            elif category in (1, 2, 3):
                draw_regions(context)
            else:
                draw_transform(context)
    finally:
        plt.close = original_close
        original_close("all")
    if not output.exists():
        raise ValueError("No output was generated")
    (destination / "request.json").write_text(json.dumps(payload, indent=2))
    (destination / "methods.txt").write_text(methods_text(context) + "\n", encoding="utf-8")
    provenance = {
        "inputs": [
            dict(
                name=r["name"],
                role=r["role"],
                group=r["group"],
                label=r["label"],
                color=r.get("color"),
                sha256=r.get("sha256") or sha256(r["path"]),
            )
            for r in records
        ],
        "software": software_versions(),
        "figure_size_inches": list(fig_size) if exact else "automatic (cropped to content)",
        "dpi": dpi,
        "synthetic_example": False,
        "example_source": "example"
        if all(r.get("source") == "example" for r in records)
        else "upload",
        "warnings": warnings,
    }
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2))
    return warnings, dict(signals=len(signals), regions=len(regions), genes=len(genes))


def draw_locus(ctx):
    options, extra = ctx["options"], ctx["extra"]
    plot, signals, regions, genes = ctx["plot"], ctx["signals"], ctx["regions"], ctx["genes"]
    track_style = extra.get("track_style") or {}
    tracks = []
    if plot != 2:
        tracks = [
            dict(
                id=name,
                type="signal",
                label=signals[index]["label"],
                source={"path": path},
                style={"color": signals[index]["color"], **track_style},
            )
            for index, (name, path) in enumerate(ctx["sources"].items())
        ]
        if plot == 1:
            tracks.append(dict(id="overlay", type="overlay", members=ctx["names"], label="Signals"))
    for region in regions:
        name = f"{region['name']} {region['label']}".lower()
        if genes and ("transcript" in name or "gene" in name):
            continue
        is_hmm = (region.get("info") or {}).get("kind") == "chromHMM" or region.get(
            "track"
        ) == "chromhmm"
        tracks.append(
            dict(
                id=f"peak_{region['id'][:8]}",
                type="chromhmm" if is_hmm else "peak",
                label=region["label"],
                source={"path": region["path"]},
                style={} if is_hmm else {"color": region["color"]},
            )
        )
    if genes:
        tracks.append(dict(id="genes", type="gene", source={"path": genes[0]["path"]}))
    fig_size = ctx["fig_size"]
    layout = {
        "width": fig_size[0] or 9,
        "height": fig_size[1],
        "track_gap": number(options, "track_gap", 0.18, 0, 0.9),
        "y_tick_count": int(number(options, "y_ticks", 3, 2, 12)),
        "y_tick_decimals": int(number(options, "y_decimals", 0, 0, 6)),
        "x_label": options.get("xlab") or None,
        "x_tick_unit": options.get("x_tick_unit", "auto"),
        "gene_model_label": options.get("gene_label", "name") or "name",
        **(extra.get("layout") or {}),
    }
    style = {
        "title": ctx["title"],
        "font_family": options.get("font_family") or "DejaVu Sans",
        "title_font": {"size": number(options, "title_size", 12, 4, 60), "weight": "bold"},
        "axis_label_font": {"size": number(options, "axis_label_size", 10, 4, 48)},
        "tick_font": {"size": number(options, "tick_size", 8, 3, 40)},
        "track_label_font": {"size": number(options, "track_label_size", 9, 3, 40)},
        "signal_fill_alpha": number(options, "fill_alpha", 0.85, 0, 1),
        "signal_linewidth": ctx["linewidth"] if options.get("linewidth") not in (None, "") else 0.6,
        **(extra.get("style") or {}),
    }
    spec = FigureSpec.parse_obj(
        dict(
            genome={"assembly": options.get("assembly") or "unknown"},
            region=dict(chrom=ctx["chrom"], start=ctx["start"], end=ctx["end"]),
            tracks=tracks,
            highlights=parse_highlights(options.get("highlights"), ctx["colors"][-1]),
            collapse_transcripts=options.get("transcripts", "collapsed") == "collapsed",
            scales={"mode": options.get("scale", "auto_highest")},
            layout=layout,
            style=style,
        )
    )
    render(spec, ctx["output"])


def _region_matrix(ctx, samples, alignment):
    options, signals, regions = ctx["options"], ctx["signals"], ctx["regions"]
    body = (
        integer(options.get("body_length", 1000) or 1000, "Body length", 10, 1_000_000)
        if alignment == "scale-regions"
        else 0
    )
    skip = boolish(options.get("skip_zeros", True))
    cache = ctx["cache"]
    key = cache.key(
        "matrix",
        signals=[_fingerprint(r["path"]) for r in signals],
        names=ctx["names"],
        regions=[_fingerprint(r["path"]) + [r["label"]] for r in regions],
        flank=ctx["flank"],
        binsize=ctx["binsize"],
        alignment=alignment,
        body=body,
        skip=skip,
    )
    cached = cache.path(key, ".npz")
    if cached and cached.is_file():
        cache.touch(cached)
        matrix = load_matrix(cached)
        ctx["notes"]["cache"] = True
    else:
        matrix = compute_region_matrix(
            samples,
            [(r["path"], r["label"]) for r in regions],
            flank=ctx["flank"],
            binsize=ctx["binsize"],
            ref_point="TSS" if alignment == "scale-regions" else alignment,
            body_length=body,
            skip_zeros=skip,
            engine="pybigwig",
        )
        if cached:
            save_matrix(matrix, cached)
    # Display metadata is not part of the cache key; refresh it from this request.
    return replace(
        matrix,
        sample_labels=tuple(s.label for s in samples),
        conditions=tuple(s.condition for s in samples),
        colormaps=tuple(s.colormap for s in samples),
        profile_colors=tuple(s.profile_color for s in samples),
    )


def _summary(ctx):
    """Sample-by-feature table over BED regions, genome bins or sliding windows (cached)."""
    from ..pipeline.genome_bins import DEFAULT_SKIP, genome_summary, region_summary

    options, cache = ctx["options"], ctx["cache"]
    mode = relationship_mode(options) if ctx["category"] == 2 else "regions"
    stat = str(options.get("rel_stat") or ("sum" if mode == "regions" else "mean"))
    skip_zeros = boolish(options.get("rel_skip_zeros", True))
    max_rows = int(number(options, "rel_max_rows", 0, 0, MAX_SUMMARY_ROWS)) or None
    params = dict(mode=mode, stat=stat, skip_zeros=skip_zeros, max_rows=max_rows)
    blacklist = next((r["path"] for r in ctx["blacklist"]), None)
    if mode == "regions":
        params["regions"] = [_fingerprint(r["path"]) for r in ctx["regions"]]
    else:
        size = integer(options.get("rel_bin") or 10000, "Bin size", 10, 100_000_000)
        step = (
            size
            if mode == "bins"
            else integer(options.get("rel_step") or size // 5, "Step", 10, size)
        )
        scope = str(options.get("rel_scope") or "genome")
        chroms = [
            c.strip()
            for c in str(options.get("rel_chroms") or "").replace(";", ",").split(",")
            if c.strip()
        ]
        skip = str(options.get("rel_skip", DEFAULT_SKIP))
        region = (ctx["chrom"], ctx["start"], ctx["end"]) if scope == "locus" else None
        params.update(
            size=size,
            step=step,
            scope=scope,
            chroms=chroms if scope == "chroms" else None,
            skip=skip,
            region=region,
            blacklist=_fingerprint(blacklist) if blacklist else None,
        )
    key = cache.key(
        "summary",
        signals=[_fingerprint(r["path"]) for r in ctx["signals"]],
        names=ctx["names"],
        **params,
    )
    cached = cache.path(key, ".pkl")
    if cached and cached.is_file():
        cache.touch(cached)
        ctx["notes"]["cache"] = True
        with cached.open("rb") as handle:
            summary, report = pickle.load(handle)  # noqa: S301 - written only by this engine
    elif mode == "regions":
        summary, report = region_summary(
            ctx["sources"],
            [r["path"] for r in ctx["regions"]],
            stat=stat,
            skip_zeros=skip_zeros,
            max_rows=max_rows,
        )
    else:
        if params["region"] is None:
            _check_window_budget(ctx["sources"], params)
        summary, report = genome_summary(
            ctx["sources"],
            size=params["size"],
            step=params["step"],
            chroms=params["chroms"],
            skip=params["skip"],
            region=params["region"],
            blacklist=blacklist,
            stat=stat,
            skip_zeros=skip_zeros,
            max_rows=max_rows,
        )
    if cached and not cached.is_file():
        with cached.open("wb") as handle:
            pickle.dump((summary, report), handle)
    if summary.values.shape[0] < 2:
        raise ValueError(
            "Fewer than two bins/regions have signal; widen the scope or switch off 'Skip all-zero rows'"
        )
    if summary.values.shape[0] > MAX_SUMMARY_ROWS:
        raise ValueError(
            f"{summary.values.shape[0]:,} rows is too many; use larger bins or 'Subsample to'"
        )
    ctx["notes"]["summary"] = dict(
        params,
        report=report.__dict__,
        blacklist_name=next((r["label"] for r in ctx["blacklist"]), None),
    )
    labels = [r["label"] for r in ctx["signals"]]
    header = "chrom\tstart\tend\t" + "\t".join(labels)
    with open(ctx["destination"] / "summary.tsv", "w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for chrom, start, end, row in zip(
            summary.chroms, summary.starts, summary.ends, summary.values, strict=True
        ):
            handle.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(f"{v:.6g}" for v in row) + "\n")
    return summary


def _check_window_budget(sources, params):
    """Refuse genome-wide requests that would allocate an unreasonable table."""
    from ..pipeline.genome_bins import SignalReader, select_chromosomes

    reader = SignalReader(next(iter(sources.values())))
    try:
        sizes = reader.sizes
    finally:
        reader.close()
    chroms = select_chromosomes(sizes, params["chroms"], params["skip"])
    windows = sum(max(0, (sizes[c] - params["size"]) // params["step"] + 1) for c in chroms)
    if windows * len(sources) > MAX_GENOME_WINDOWS * 4:
        raise ValueError(
            f"{windows:,} windows × {len(sources)} samples is too large for one figure; use larger bins, "
            "a larger step, or restrict chromosomes"
        )


def draw_regions(ctx):
    category, plot, options, extra = ctx["category"], ctx["plot"], ctx["options"], ctx["extra"]
    signals, regions, names, colors = ctx["signals"], ctx["regions"], ctx["names"], ctx["colors"]
    destination, output = ctx["destination"], ctx["output"]
    table_only = category == 2 or (category == 3 and plot == 0)
    if not (category == 2 and relationship_mode(options) != "regions"):
        region_rows = sum(len(read_bed(r["path"])) for r in regions)
        limit = MAX_SUMMARY_ROWS if table_only else MAX_REGIONS
        if not 2 <= region_rows <= limit:
            raise ValueError(
                f"Choose between 2 and {limit:,} regions per figure (got {region_rows:,})"
            )
        ctx["notes"]["regions"] = region_rows
    labels = tuple(r["label"] for r in signals)
    # Show plain labels; fall back to the internal unique names only when labels repeat.
    display_names = tuple(labels) if len(set(labels)) == len(labels) else tuple(names)
    if table_only:
        summary = _summary(ctx)
        if category == 2:
            if plot == 0:
                table = (
                    np.log1p(np.clip(summary.values, 0, None))
                    if options.get("corr_transform") == "log1p"
                    else summary.values
                )
                correlation_heatmap(
                    table,
                    labels,
                    output,
                    method=options.get("method", "pearson"),
                    cmap=resolve_cmap(options.get("corr_cmap"), "RdBu_r"),
                    annotate=boolish(options.get("annotate", True)),
                )
            elif plot == 1:
                draw_pca(ctx, summary, display_names)
            else:
                pair_grid(
                    summary.values,
                    labels,
                    output,
                    max_rows=integer(options.get("maxrows", 1000), "Rows", 10, 50_000),
                    color=colors[0],
                    hist_color=tint(colors[0], 0.6),
                )
            return
        conditions = {name: r["group"] for name, r in zip(names, signals, strict=True)}
        numerator = str(options.get("numerator") or "Treatment")
        denominator = str(options.get("denominator") or "Control")
        counts = Counter(conditions.values())
        if numerator == denominator or counts[numerator] < 2 or counts[denominator] < 2:
            raise ValueError(
                "Volcano plots need two distinct groups with ≥2 replicates each. "
                f"Current signal groups: {dict(counts)}. Set each signal's group, then pick "
                "numerator/denominator."
            )
        result = diffpeak(summary, conditions, numerator=numerator, denominator=denominator)
        volcano_plot(
            result,
            output=output,
            fdr=number(options, "fdr", 0.05, 0, 1),
            up_color=valid_color(options.get("up_color"), colors[1 % len(colors)]),
            down_color=valid_color(options.get("down_color"), colors[0]),
        )
        np.savetxt(
            destination / "statistics.tsv",
            np.column_stack([result.log2_fold_change, result.p_value, result.fdr]),
            delimiter="\t",
            header="log2_fold_change\tp_value\tfdr",
        )
        ctx["warnings"].append(
            "diffpeak uses a normal-tail approximation for its test statistic; use validated "
            "differential-analysis methods (DESeq2, edgeR, DiffBind) for inferential reporting."
        )
        return
    samples = [
        SampleMeta(
            path,
            name,
            r["group"],
            "signals",
            cmap_spec(r.get("cmap") or options.get("cmap"), r["color"]),
            r["color"],
        )
        for (name, path), r in zip(ctx["sources"].items(), signals, strict=True)
    ]
    alignment = options.get("alignment", "center")
    if alignment not in {"center", "TSS", "TES", "scale-regions"}:
        raise ValueError("Alignment must be center, TSS, TES or scale-regions")
    matrix = _region_matrix(ctx, samples, alignment)
    if category == 3:
        peak_center_comparisons(
            matrix,
            destination / "comparison",
            flank_bins=integer(options.get("center_bins", 5), "Center bins", 1, 200),
            cmap=resolve_cmap(options.get("pc_cmap"), "viridis_r"),
        )
        ctx["warnings"].append("Peak-center p-values are unadjusted Mann–Whitney tests.")
        return
    if options.get("normalization", "raw") not in ("raw", "", None):
        matrix = renormalize(matrix, approach=options["normalization"])
    sort_by = (
        options.get("sort_by")
        if options.get("sort_by") in {"mean", "median", "max", "min", "sum"}
        else "mean"
    )
    if options.get("sort", "descend") != "none":
        matrix = sort_matrix(
            matrix, method=sort_by, descending=options.get("sort", "descend") == "descend"
        )
    cluster_method = options.get("cluster", "none")
    if plot == 3 and cluster_method in {"none", "", None}:
        cluster_method = "kmeans"
    if cluster_method not in {"none", "", None}:
        k = integer(options.get("clusters", 3), "Clusters", 2, min(20, len(matrix.coords)))
        labels_ = cluster_labels(matrix, method=cluster_method, n_clusters=k, seed=7201)
        order = np.argsort(labels_, kind="stable")
        unique, counts_ = np.unique(labels_, return_counts=True)
        matrix = replace(
            matrix,
            values=matrix.values[order],
            coords=[matrix.coords[i] for i in order],
            group_labels=tuple(f"C{i + 1}" for i in unique),
            group_boundaries=tuple([0, *np.cumsum(counts_).tolist()]),
        )
        ctx["notes"]["cluster"] = (cluster_method, k)
    save_matrix(matrix, destination / "matrix.npz")
    heatmap_style, profile_style = plot_styles(
        options, ctx["title"], show_profile=plot == 2, extra=extra
    )
    heatmap_style.labels = labels
    profile_style.labels = labels
    profile_style.colors = tuple(r["color"] for r in signals)
    if ctx["fig_size"][0]:
        heatmap_style.width_per_sample = max(0.8, ctx["fig_size"][0] / max(1, len(signals)))
    if plot == 1:
        plot_profile_lines(matrix, output, style=profile_style)
    else:
        plot_enrichment_heatmap(matrix, output, style=heatmap_style)


MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
PLOTLY_SYMBOLS = [
    "circle",
    "square",
    "triangle-up",
    "diamond",
    "triangle-down",
    "cross",
    "x",
    "star",
]


def repel_labels(axis, xs, ys, texts, *, fontsize, color="#1f2933", iterations=300):
    """Place point labels without overlaps (a small ggrepel). Returns offsets in points.

    Labels start up-right of their point and are pushed apart from each other and from all
    markers in display space; labels that travel far get a thin leader line.
    """
    figure = axis.figure
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    scale = figure.dpi / 72.0
    points = axis.transData.transform(np.column_stack([xs, ys]))
    offsets = np.tile([6.0, 6.0], (len(texts), 1))  # points
    notes = [
        axis.annotate(
            t,
            (x, y),
            xytext=tuple(o),
            textcoords="offset points",
            fontsize=fontsize,
            color=color,
            ha="left",
            va="bottom",
            zorder=4,
        )
        for t, x, y, o in zip(texts, xs, ys, offsets, strict=True)
    ]
    sizes = (
        np.array([[b.width, b.height] for b in (n.get_window_extent(renderer) for n in notes)]) + 4
    )
    box = axis.get_window_extent(renderer)
    radius = 7 * scale
    for _ in range(iterations):
        lower = points + offsets * scale  # lower-left corner of each label (pixels)
        centers = lower + sizes / 2
        push = np.zeros_like(offsets)
        moved = False
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                overlap = (sizes[i] + sizes[j]) / 2 - np.abs(centers[i] - centers[j])
                if (overlap > 0).all():
                    direction = centers[i] - centers[j]
                    axis_ = int(np.argmin(overlap))  # separate along the cheaper axis
                    sign = 1.0 if direction[axis_] >= 0 else -1.0
                    step = np.zeros(2)
                    step[axis_] = sign * (overlap[axis_] / 2 + 1)
                    push[i] += step
                    push[j] -= step
                    moved = True
            for k in range(len(points)):  # keep labels off every marker
                near = np.clip(points[k], lower[i], lower[i] + sizes[i])
                gap = points[k] - near
                distance = float(np.hypot(*gap))
                if distance < radius:
                    away = (centers[i] - points[k]) / (np.hypot(*(centers[i] - points[k])) or 1)
                    push[i] += away * (radius - distance + 1)
                    moved = True
        offsets += push / scale * 0.5
        lower = points + offsets * scale
        offsets[:, 0] += np.clip(box.x0 - lower[:, 0], 0, None) / scale  # stay inside the axes
        offsets[:, 0] -= np.clip(lower[:, 0] + sizes[:, 0] - box.x1, 0, None) / scale
        offsets[:, 1] += np.clip(box.y0 - lower[:, 1], 0, None) / scale
        offsets[:, 1] -= np.clip(lower[:, 1] + sizes[:, 1] - box.y1, 0, None) / scale
        if not moved:
            break
    for note in notes:
        note.remove()
    for text, x, y, offset in zip(texts, xs, ys, offsets, strict=True):
        far = np.hypot(*offset) > 16  # leader line only when the label had to move away
        axis.annotate(
            text,
            (x, y),
            xytext=tuple(offset),
            textcoords="offset points",
            fontsize=fontsize,
            color=color,
            ha="left",
            va="bottom",
            zorder=4,
            arrowprops=dict(arrowstyle="-", color="#9aa5b1", lw=0.6, shrinkA=0, shrinkB=4)
            if far
            else None,
        )
    return offsets


def draw_pca(ctx, summary, display_names):
    """PCA scatter with repelled labels, group markers, and an interactive companion (pca.json)."""
    options, signals, destination = ctx["options"], ctx["signals"], ctx["destination"]
    top = integer(options.get("pca_top") or 1000, "Top features", 2, MAX_SUMMARY_ROWS)
    log2 = options.get("log", "log2") == "log2"
    result = pca_scores(replace(summary, samples=display_names), log2=log2, top=top, n_components=3)
    coords = np.asarray(result.coordinates, dtype=float)
    variance = np.asarray(result.variance_explained, dtype=float)
    pct = [100 * variance[i] if i < variance.size else 0.0 for i in range(3)]
    groups = [r["group"] or "" for r in signals]
    group_order = list(dict.fromkeys(groups))
    use_groups = 1 < len(group_order) < len(groups)  # only when groups actually pool samples
    size = number(options, "point_size", 60, 4, 800)
    fig_size = ctx["fig_size"]
    figure, axis = plt.subplots(figsize=(fig_size[0] or 5.6, fig_size[1] or 4.6))
    for i in range(len(display_names)):
        marker = MARKERS[group_order.index(groups[i]) % len(MARKERS)] if use_groups else "o"
        axis.scatter(
            coords[i, 0],
            coords[i, 1],
            s=size,
            color=signals[i]["color"],
            marker=marker,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
    axis.axhline(0, color="#e3e7ec", linewidth=0.8, zorder=1)
    axis.axvline(0, color="#e3e7ec", linewidth=0.8, zorder=1)
    axis.set_xlabel(f"PC1 ({pct[0]:.1f}%)", fontsize=number(options, "axis_label_size", 10, 4, 48))
    axis.set_ylabel(f"PC2 ({pct[1]:.1f}%)", fontsize=number(options, "axis_label_size", 10, 4, 48))
    axis.tick_params(labelsize=number(options, "tick_size", 8, 3, 40))
    axis.margins(x=0.22, y=0.18)
    if use_groups:
        from matplotlib.lines import Line2D

        handles = [
            Line2D(
                [],
                [],
                marker=MARKERS[k % len(MARKERS)],
                linestyle="",
                color="#52606d",
                label=g or "(none)",
            )
            for k, g in enumerate(group_order)
        ]
        axis.legend(
            handles=handles,
            frameon=False,
            fontsize=number(options, "legend_size", 9, 3, 40),
            loc="best",
            title="group",
            title_fontsize=number(options, "legend_size", 9, 3, 40),
        )
    axis.set_title(ctx["title"], fontsize=number(options, "title_size", 12, 4, 60), pad=10)
    figure.tight_layout()
    offsets = repel_labels(
        axis,
        coords[:, 0],
        coords[:, 1],
        list(display_names),
        fontsize=number(options, "sample_label_size", 9, 3, 40),
    )
    data = dict(
        title=ctx["title"],
        variance=pct,
        top=top,
        log2=log2,
        font=options.get("font_family") or "",
        points=[
            dict(
                name=name,
                group=groups[i],
                color=signals[i]["color"],
                file=signals[i]["name"],
                pcs=[float(coords[i, k]) if k < coords.shape[1] else 0.0 for k in range(3)],
                symbol=PLOTLY_SYMBOLS[group_order.index(groups[i]) % len(PLOTLY_SYMBOLS)]
                if use_groups
                else "circle",
                label_offset=[float(offsets[i, 0]), float(offsets[i, 1])],
            )
            for i, name in enumerate(display_names)
        ],
        point_size=size,
        use_groups=use_groups,
    )
    (destination / "pca.json").write_text(json.dumps(data, indent=1))
    with open(destination / "pca_scores.tsv", "w", encoding="utf-8") as handle:
        handle.write("sample\tgroup\tfile\tPC1\tPC2\tPC3\n")
        for point in data["points"]:
            handle.write(
                "\t".join(
                    [
                        point["name"],
                        point["group"],
                        point["file"],
                        *(f"{v:.6g}" for v in point["pcs"]),
                    ]
                )
                + "\n"
            )
    _pca_html(data, destination / "pca_interactive.html")
    plt.close(figure)


def _pca_html(data, path):
    """Standalone interactive PCA (Plotly from its CDN) for sharing outside the studio."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        return
    v = data["variance"]
    fig = go.Figure()
    for point in data["points"]:
        fig.add_trace(
            go.Scatter(
                x=[point["pcs"][0]],
                y=[point["pcs"][1]],
                mode="markers",
                name=point["name"],
                marker=dict(
                    size=13,
                    color=point["color"],
                    symbol=point["symbol"],
                    line=dict(color="white", width=1),
                ),
                customdata=[[point["group"], point["pcs"][2], point["file"]]],
                hovertemplate=(
                    f"<b>{point['name']}</b><br>group: %{{customdata[0]}}<br>PC1: %{{x:.2f}} ({v[0]:.1f}%)"
                    f"<br>PC2: %{{y:.2f}} ({v[1]:.1f}%)<br>PC3: %{{customdata[1]:.2f}} ({v[2]:.1f}%)"
                    "<br>%{customdata[2]}<extra></extra>"
                ),
            )
        )
        fig.add_annotation(
            x=point["pcs"][0],
            y=point["pcs"][1],
            text=point["name"],
            showarrow=True,
            arrowcolor="#9aa5b1",
            arrowwidth=0.8,
            ax=point["label_offset"][0] * 1.33 + 10,
            ay=-point["label_offset"][1] * 1.33 - 6,
            font=dict(size=12, color="#1f2933"),
        )
    fig.update_layout(
        title=data["title"],
        template="simple_white",
        showlegend=False,
        hovermode="closest",
        xaxis_title=f"PC1 ({v[0]:.1f}%)",
        yaxis_title=f"PC2 ({v[1]:.1f}%)",
        xaxis=dict(zeroline=True, zerolinecolor="#e3e7ec"),
        yaxis=dict(zeroline=True, zerolinecolor="#e3e7ec"),
    )
    fig.write_html(
        str(path),
        include_plotlyjs="cdn",
        full_html=True,
        config=dict(
            editable=False,
            edits=dict(annotationTail=True),
            displaylogo=False,
            toImageButtonOptions=dict(format="svg", filename="pca"),
        ),
    )


def draw_transform(ctx):
    plot, options, signals = ctx["plot"], ctx["options"], ctx["signals"]
    chrom, start, end, flank = ctx["chrom"], ctx["start"], ctx["end"], ctx["flank"]
    pixels = integer(options.get("pixels", 800) or 800, "Resolution", 50, 5000)
    a_rec, b_rec = signals[:2]
    if options.get("direction", "b-a") == "a-b":
        a_rec, b_rec = b_rec, a_rec
    a, b = (
        query_signal(r["path"], chrom, start, end, pixels=pixels, engine="pybigwig")
        for r in (a_rec, b_rec)
    )
    pseudo = float(options.get("pseudocount", 0.01) or 0.01)
    if not np.isfinite(pseudo) or pseudo <= 0:
        raise ValueError("Pseudocount must be positive")
    name = CATALOG[4][1][plot]
    if plot == 3:
        path = b_rec["path"]
        left = (
            query_signal(
                path, chrom, max(0, start - flank), start, pixels=100, engine="pybigwig"
            ).values
            if start
            else np.array([])
        )
        right = query_signal(path, chrom, end, end + flank, pixels=100, engine="pybigwig").values
        background = np.concatenate([left, right])
        background = background[np.isfinite(background)]
        if not background.size:
            raise ValueError(
                "No finite background data in flanking windows; increase Flank or move the window"
            )
        estimator = options.get("estimator", "median")
        baseline = float(np.median(background) if estimator == "median" else np.mean(background))
        values = transform_values(a.values, baseline, "subtract")
    else:
        values = compare_signal(
            a, b, operation=["difference", "ratio", "log2_ratio"][plot], pseudocount=pseudo
        ).values
    derived = valid_color(options.get("derived_color"), "#3b4252")
    fig_size = ctx["fig_size"]
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(fig_size[0] or 9, fig_size[1] or 6),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1, 1.25]},
    )
    x = a.bin_starts
    for ax, vals, label, col in zip(
        axes,
        [a.values, b.values, values],
        [a_rec["label"], b_rec["label"], name],
        [a_rec["color"], b_rec["color"], derived],
        strict=True,
    ):
        ax.plot(x, vals, color=col, linewidth=ctx["linewidth"])
        baseline_value = 1.0 if plot == 1 else 0.0
        ax.fill_between(
            x,
            baseline_value if ax is axes[-1] else 0,
            vals,
            color=col,
            alpha=ctx["fill_alpha"],
            linewidth=0,
        )
        if ax is axes[-1]:
            ax.axhline(baseline_value, color="#9aa5b1", linewidth=0.7, linestyle="--")
        ax.set_ylabel(label, fontsize=number(options, "axis_label_size", 10, 4, 48))
        ax.tick_params(labelsize=number(options, "tick_size", 8, 3, 40))
    axes[-1].set_xlabel(options.get("xlab") or f"{chrom} coordinate (bp)")
    np.savetxt(
        ctx["destination"] / "signals.tsv",
        np.column_stack([a.bin_starts, a.values, b.values, values]),
        delimiter="\t",
        header="start\tA\tB\tderived",
    )
    plt.close(fig)


# ---------------------------------------------------------------- methods text


def _bp(value):
    value = int(value)
    return f"{value / 1000:g} kb" if value >= 1000 and value % 100 == 0 else f"{value} bp"


def methods_text(ctx) -> str:
    """Draft Methods paragraph describing exactly how this figure was produced."""
    category, plot, options = ctx["category"], ctx["plot"], ctx["options"]
    versions = software_versions()
    gc = versions.get("genomecanvas", "unversioned")
    mpl = versions.get("matplotlib", "?")
    n_sig = len(ctx["signals"])
    sig = f"{n_sig} signal track{'s' if n_sig != 1 else ''} (BigWig/bedGraph)"
    parts = []
    if category == 0:
        parts.append(
            f"Genome browser-style tracks for {ctx['chrom']}:{ctx['start'] + 1:,}-{ctx['end']:,} "
            f"({options.get('assembly') or 'assembly not specified'}) were drawn from {sig}"
            + (", with gene models from a GTF annotation" if ctx["genes"] else "")
            + ". Signal was summarised per display bin; y-axes were "
            + (
                "shared across tracks."
                if options.get("scale", "auto_highest") == "auto_highest"
                else "scaled per track."
            )
        )
    elif category in (1, 3) and not (category == 3 and plot == 0):
        alignment = options.get("alignment", "center")
        anchor = {
            "center": "region centres",
            "TSS": "region starts (TSS for stranded features)",
            "TES": "region ends (TES)",
            "scale-regions": "scaled region bodies",
        }[alignment]
        text = (
            f"Signal from {sig} was extracted around {anchor} of {ctx['notes'].get('regions', 'the')} "
            f"regions ({len(ctx['regions'])} region set{'s' if len(ctx['regions']) != 1 else ''}) "
            f"in ±{_bp(ctx['flank'])} windows using {_bp(ctx['binsize'])} bins"
        )
        if alignment == "scale-regions":
            text += f", with region bodies scaled to {_bp(options.get('body_length', 1000))}"
        if boolish(options.get("skip_zeros", True)):
            text += "; regions with no signal were excluded"
        text += "."
        parts.append(text)
        if category == 1:
            norm = options.get("normalization", "raw")
            if norm not in ("raw", "", None):
                parts.append(
                    "Matrices were "
                    + (
                        "min–max scaled per sample."
                        if norm == "profileMinMax"
                        else "quantile normalised."
                    )
                )
            if options.get("plot_scale", "linear") != "linear":
                parts.append(f"Values are displayed on a {options['plot_scale']} scale.")
            if "cluster" in ctx["notes"]:
                method, k = ctx["notes"]["cluster"]
                parts.append(
                    f"Rows were grouped into {k} clusters by {'k-means' if 'kmeans' in method else 'hierarchical clustering'} (seed 7201)."
                )
            elif options.get("sort", "descend") != "none":
                parts.append(
                    f"Rows were sorted by {options.get('sort_by', 'mean')} signal ({options.get('sort', 'descend')}ing)."
                )
            if plot in (1, 2) and boolish(options.get("error_estimates", True)):
                parts.append("Profiles show the mean ± standard error of the mean.")
        else:
            parts.append(
                f"Mean signal within ±{options.get('center_bins', 5)} bins of the centre was compared between "
                "region sets with two-sided Mann–Whitney U tests (p-values not adjusted)."
            )
    elif category == 2:
        parts.append(_summary_methods(ctx, sig))
        if plot == 0:
            parts.append(
                f"{options.get('method', 'pearson').capitalize()} correlation coefficients between samples were computed"
                + (
                    " on log(x + 1)-transformed values"
                    if options.get("corr_transform") == "log1p"
                    else ""
                )
                + " and samples ordered by hierarchical clustering (Ward)."
            )
        elif plot == 1:
            parts.append(
                "Principal component analysis was performed on "
                + ("log2-transformed " if options.get("log", "log2") == "log2" else "")
                + "values of the most variable features."
            )
        else:
            parts.append(
                "Pairwise scatter plots show log10(signal + 1); upper panels give Pearson r."
            )
    elif category == 3:
        parts.append(
            f"Signal from {sig} was summed over each region and compared between "
            f"{options.get('numerator')} and {options.get('denominator')} using GenomeCanvas diffpeak "
            f"(normal approximation; FDR < {options.get('fdr', 0.05)}). This is exploratory; formal "
            "differential analysis should use dedicated count-based methods."
        )
    else:
        op = [
            "difference (A − B)",
            "ratio (A / B)",
            "log2 ratio (log2[(A + p)/(B + p)])",
            "background-subtracted signal",
        ][plot]
        parts.append(
            f"The {op} was computed for {ctx['chrom']}:{ctx['start'] + 1:,}-{ctx['end']:,} at "
            f"{options.get('pixels', 800)} bins"
            + (
                f" with pseudocount p = {options.get('pseudocount', 0.01)}"
                if plot in (1, 2)
                else ""
            )
            + (
                f"; background was the {options.get('estimator', 'median')} of ±{_bp(ctx['flank'])} flanking windows"
                if plot == 3
                else ""
            )
            + "."
        )
    parts.append(f"Figures were generated with GenomeCanvas {gc} (Matplotlib {mpl}).")
    return " ".join(parts)


def _summary_methods(ctx, sig):
    info = ctx["notes"].get("summary") or {}
    report = info.get("report") or {}
    stat = {
        "mean": "mean signal over covered bases",
        "mean0": "mean signal (uncovered bases = 0)",
        "sum": "summed signal (signal × bp)",
        "max": "maximum signal",
    }.get(info.get("stat"), "signal")
    if info.get("mode") in ("bins", "windows"):
        if info["mode"] == "bins":
            unit = f"non-overlapping {_bp(info['size'])} genomic bins"
        else:
            unit = f"sliding {_bp(info['size'])} windows with a {_bp(info['step'])} step"
        where = (
            f"within {info['region'][0]}:{info['region'][1] + 1:,}-{info['region'][2]:,}"
            if info.get("region")
            else f"across {len(report.get('chroms', []))} chromosomes"
            + (f" (excluding {info['skip']})" if info.get("skip") else "")
        )
        text = f"The {stat} of {sig} was computed in {unit} {where}"
        if info.get("blacklist_name"):
            text += f"; windows overlapping {info['blacklist_name']} were removed ({report.get('blacklisted', 0):,})"
    else:
        text = f"The {stat} of {sig} was computed for each of {report.get('windows_total', 'the'):,} regions"
    text += "."
    if info.get("skip_zeros"):
        text += f" Features with zero signal in every sample were excluded ({report.get('windows_kept', 0):,} retained)."
    if report.get("subsampled"):
        text += f" A random subset of {info.get('max_rows'):,} features (seed 7201) was used."
    return text


def meta(extra: dict[str, Any] | None = None) -> dict:
    """Static catalog/style metadata for the front end."""
    import matplotlib.font_manager as font_manager

    available = {f.name for f in font_manager.fontManager.ttflist}

    def stops(name):
        cmap = resolve_cmap(name if name != "auto" else "#ffffff,#2A7F78,#14403c", "Blues")
        return [to_hex(cmap(i / 6)) for i in range(7)]

    return dict(
        engine="GenomeCanvas",
        catalog=CATALOG,
        requirements=REQUIREMENTS,
        palettes={
            k: dict(name=v[0], colors=v[1])
            for k, v in PALETTES.items()
            if k not in ("teal", "blue")
        },
        colormaps={name: stops(name) for name in COLORMAPS},
        diverging={name: stops(name) for name in DIVERGING},
        themes=list(THEMES),
        fonts=[f for f in FONTS if f in available or f == "DejaVu Sans"],
        suffixes=sorted(SUFFIXES),
        style_fields=dict(
            heatmap=[f.name for f in dataclasses.fields(HeatmapPlotStyle)],
            profile=[f.name for f in dataclasses.fields(ProfilePlotStyle)],
        ),
        **(extra or {}),
    )
