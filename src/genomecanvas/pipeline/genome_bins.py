"""Fast sample-by-feature signal tables: genome-wide bins, sliding windows and BED regions.

Equivalent in purpose to deepTools ``multiBigwigSummary bins`` / ``BED-file``. Two read paths
are chosen per file by data density (BigWig bytes per genomic base):

* sparse/run-length data (bamCoverage bins, broad zero runs): raw intervals are read per
  chromosome chunk and integrated with cumulative sums, so cost scales with the number of
  intervals, not genome length;
* dense data (bp-resolution signal): per-base values are read as numpy arrays in small chunks
  and reduced with reshape, keeping memory bounded and avoiding millions of Python tuples.

Statistics (per bin, window or region):

* ``mean``  – mean over bases *with data* (pyBigWig/deepTools convention); 0 if no data.
* ``mean0`` – mean over all bases, missing bases counted as 0.
* ``sum``   – bp-weighted sum (signal × bases), missing = 0. Matches ``extract_summary``.
* ``max``   – maximum value.
"""

from __future__ import annotations

import fnmatch
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..analysis.summary import RegionSummary
from ..io.bed import read_bed

DEFAULT_SKIP = "chrM, MT, *_random, chrUn*, *_alt, *_fix, chrEBV, *hap*, *decoy*"
STATS = ("mean", "mean0", "sum", "max")
CHUNK_BP = 20_000_000  # interval path
VALUES_CHUNK_BP = 4_000_000  # per-base path: ~32 MB of float64 per chunk
DENSE_BYTES_PER_BASE = 0.1  # 1-bp BigWigs are ~5 B/base; 10 bp ~0.5; sparse bins ~0.005


class SignalReader:
    """Interval access for BigWig or bedGraph files."""

    def __init__(self, path: str | Path, read_mode: str = "auto"):
        self.path = str(path)
        self.bigwig = self.path.lower().endswith((".bw", ".bigwig"))
        self.dense = False
        if self.bigwig:
            import pyBigWig

            self.handle = pyBigWig.open(self.path)
            if self.handle is None:
                raise ValueError(f"Cannot open BigWig {self.path}")
            self.sizes = dict(self.handle.chroms())
            density = Path(self.path).stat().st_size / max(1, sum(self.sizes.values()))
            self.dense = read_mode == "values" or (
                read_mode == "auto" and density > DENSE_BYTES_PER_BASE
            )
            self.dense = self.dense and bool(pyBigWig.numpy)
        else:
            self._load_bedgraph()

    def _load_bedgraph(self):
        rows: dict[str, list] = {}
        with open(self.path, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                chrom, start, end, value = line.split()[:4]
                rows.setdefault(chrom, []).append((int(start), int(end), float(value)))
        self.data = {}
        for chrom, items in rows.items():
            array = np.array(sorted(items), dtype=float)
            self.data[chrom] = (array[:, 0], array[:, 1], array[:, 2])
        self.sizes = {chrom: int(arrays[1].max()) for chrom, arrays in self.data.items()}

    def intervals(self, chrom: str, start: int, end: int):
        """Sorted (starts, ends, values) clipped to [start, end); NaN values dropped."""
        if self.bigwig:
            raw = (
                self.handle.intervals(chrom, int(start), int(end)) if chrom in self.sizes else None
            )
            if not raw:
                return _empty()
            array = np.asarray(raw, dtype=float)
            s, e, v = array[:, 0], array[:, 1], array[:, 2]
        else:
            if chrom not in self.data:
                return _empty()
            s, e, v = self.data[chrom]
            lo = np.searchsorted(e, start, side="right")
            hi = np.searchsorted(s, end, side="left")
            s, e, v = s[lo:hi], e[lo:hi], v[lo:hi]
        keep = np.isfinite(v)
        s, e, v = np.clip(s[keep], start, end), np.clip(e[keep], start, end), v[keep]
        nonempty = e > s
        return s[nonempty], e[nonempty], v[nonempty]

    def values(self, chrom: str, start: int, end: int):
        """Per-base values (NaN = no data) for [start, end) of a BigWig."""
        start, end = int(start), int(end)
        out = np.full(end - start, np.nan)
        length = self.sizes.get(chrom)
        if length is None:
            return out
        lo, hi = max(0, start), min(end, length)  # windows may run past a chromosome edge
        if hi > lo:
            out[lo - start : hi - start] = self.handle.values(chrom, lo, hi, numpy=True)
        return out

    def close(self):
        if self.bigwig:
            self.handle.close()


def _empty():
    z = np.zeros(0)
    return z, z, z


def _integrals(s, e, v, positions):
    """∫signal and ∫coverage from the first interval to each position (piecewise linear)."""
    if not s.size:
        zero = np.zeros(positions.shape)
        return zero, zero
    width = e - s
    cum_sum = np.concatenate([[0.0], np.cumsum(v * width)])
    cum_cov = np.concatenate([[0.0], np.cumsum(width)])
    k = np.searchsorted(s, positions, side="right") - 1
    inside = k >= 0
    kk = np.where(inside, k, 0)
    part = np.clip(positions - s[kk], 0, width[kk])
    total = np.where(inside, cum_sum[kk] + v[kk] * part, 0.0)
    cover = np.where(inside, cum_cov[kk] + part, 0.0)
    return total, cover


def _bin_max(s, e, v, lo, step, n):
    out = np.full(n, -np.inf)
    if s.size:
        b0 = ((s - lo) // step).astype(np.int64)
        b1 = ((e - 1 - lo) // step).astype(np.int64)
        span = b1 - b0 + 1
        bins = np.repeat(b0, span) + (
            np.arange(span.sum()) - np.repeat(np.cumsum(span) - span, span)
        )
        values = np.repeat(v, span)
        ok = (bins >= 0) & (bins < n)
        np.maximum.at(out, bins[ok], values[ok])
    return out


def bin_chromosome(reader: SignalReader, chrom: str, lo: int, hi: int, step: int, stat: str):
    """Per-step-bin (sum, coverage, max) arrays over [lo, lo + n*step)."""
    n = (hi - lo) // step
    total = np.zeros(n)
    cover = np.zeros(n)
    peak = np.full(n, -np.inf) if stat == "max" else None
    if reader.dense:
        chunk = max(step, (VALUES_CHUNK_BP // step) * step)
        for first in range(0, n, chunk // step):
            last = min(n, first + chunk // step)
            block = reader.values(chrom, lo + first * step, lo + last * step).reshape(
                last - first, step
            )
            finite = np.isfinite(block)
            total[first:last] = np.where(finite, block, 0.0).sum(axis=1)
            cover[first:last] = finite.sum(axis=1)
            if peak is not None:
                peak[first:last] = np.where(finite, block, -np.inf).max(axis=1)
        return total, cover, peak
    chunk = max(step, (CHUNK_BP // step) * step)
    for first in range(0, n, chunk // step):
        last = min(n, first + chunk // step)
        a, b = lo + first * step, lo + last * step
        s, e, v = reader.intervals(chrom, a, b)
        edges = a + step * np.arange(last - first + 1, dtype=float)
        f_sum, f_cov = _integrals(s, e, v, edges)
        total[first:last] = np.diff(f_sum)
        cover[first:last] = np.diff(f_cov)
        if peak is not None:
            peak[first:last] = _bin_max(s, e, v, a, step, last - first)
    return total, cover, peak


def _finish(stat, total, cover, peak, width):
    if stat == "sum":
        return total
    if stat == "mean0":
        return total / width
    if stat == "max":
        return np.where(np.isfinite(peak), peak, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cover > 0, total / np.maximum(cover, 1e-12), 0.0)


def _rolling(values, w, how):
    if w == 1:
        return values
    if how == "max":
        from numpy.lib.stride_tricks import sliding_window_view

        return sliding_window_view(values, w).max(axis=1) if values.size >= w else values[:0]
    c = np.concatenate([[0.0], np.cumsum(values)])
    return c[w:] - c[:-w]


def select_chromosomes(sizes: dict[str, int], include=None, skip=DEFAULT_SKIP) -> list[str]:
    def norm(name):
        return name.lower().removeprefix("chr")

    chosen = list(sizes)
    if include:
        wanted = {norm(c) for c in include}
        chosen = [c for c in chosen if norm(c) in wanted]
    patterns = [p.strip().lower() for p in str(skip or "").split(",") if p.strip()]
    chosen = [c for c in chosen if not any(fnmatch.fnmatch(c.lower(), p) for p in patterns)]

    def order(name):  # natural order: chr1, chr2, … chr10, chrX
        key = norm(name)
        return (0, int(key), "") if key.isdigit() else (1, 0, key)

    return sorted(chosen, key=order)


@dataclass
class BinReport:
    windows_total: int = 0
    windows_kept: int = 0
    chroms: list[str] = field(default_factory=list)
    blacklisted: int = 0
    zero_removed: int = 0
    subsampled: bool = False
    seconds: float = 0.0


def genome_summary(
    sources: dict[str, str],
    *,
    size: int = 10_000,
    step: int | None = None,
    chroms: list[str] | None = None,
    skip: str = DEFAULT_SKIP,
    region: tuple[str, int, int] | None = None,
    blacklist: str | Path | None = None,
    stat: str = "mean",
    skip_zeros: bool = True,
    max_rows: int | None = None,
    seed: int = 7201,
) -> tuple[RegionSummary, BinReport]:
    """Signal table over fixed bins (step = size) or sliding windows (step < size)."""
    began = time.time()
    if stat not in STATS:
        raise ValueError(f"Statistic must be one of {', '.join(STATS)}")
    step = int(step or size)
    if size < 10 or step < 10:
        raise ValueError("Bin size and step must be at least 10 bp")
    if step > size:
        raise ValueError("Step cannot be larger than the window size")
    if size % step:
        raise ValueError(
            "Window size must be a multiple of the step (e.g. 10 kb windows, 2 kb step)"
        )
    w = size // step
    readers = {label: SignalReader(path) for label, path in sources.items()}
    try:
        shared = None
        for reader in readers.values():
            shared = (
                dict(reader.sizes)
                if shared is None
                else {c: min(n, reader.sizes[c]) for c, n in shared.items() if c in reader.sizes}
            )
        if not shared:
            raise ValueError("The selected signals share no chromosome names")
        if region:
            chrom, lo, hi = region
            match = [c for c in shared if c == chrom] or [
                c
                for c in shared
                if c.lower().removeprefix("chr") == chrom.lower().removeprefix("chr")
            ]
            if not match:
                raise ValueError(f"{chrom} is not present in all selected signals")
            spans = {match[0]: (max(0, int(lo)), min(int(hi), shared[match[0]]))}
        else:
            spans = {c: (0, shared[c]) for c in select_chromosomes(shared, chroms, skip)}
        if not spans:
            raise ValueError("No chromosomes left after applying the include/skip lists")
        black = {}
        if blacklist:
            for feature in read_bed(blacklist):
                black.setdefault(feature.interval.chrom, []).append(
                    (feature.interval.start, feature.interval.end)
                )
        report = BinReport(chroms=list(spans))
        columns: list[list[np.ndarray]] = [[] for _ in readers]
        coords_chrom, coords_start = [], []
        for chrom, (lo, hi) in spans.items():
            n_steps = (hi - lo) // step
            n_win = n_steps - w + 1
            if n_win <= 0:
                continue
            for index, reader in enumerate(readers.values()):
                total, cover, peak = bin_chromosome(
                    reader, chrom, lo, lo + n_steps * step, step, stat
                )
                if stat == "max":
                    values = _rolling(np.where(np.isfinite(peak), peak, 0.0), w, "max")
                else:
                    t, c = _rolling(total, w, "sum"), _rolling(cover, w, "sum")
                    values = _finish(stat, t, c, None, size)
                columns[index].append(values[:n_win])
            starts = lo + step * np.arange(n_win, dtype=np.int64)
            keep = np.ones(n_win, dtype=bool)
            for b_start, b_end in black.get(chrom, []):
                k_min = max(0, int(np.floor((b_start - size - lo) / step)) + 1)
                k_max = min(n_win - 1, int(np.ceil((b_end - lo) / step)) - 1)
                if k_max >= k_min:
                    keep[k_min : k_max + 1] = False
            report.blacklisted += int((~keep).sum())
            for index in range(len(readers)):
                columns[index][-1] = columns[index][-1][keep]
            coords_chrom.append(np.full(int(keep.sum()), chrom, dtype=object))
            coords_start.append(starts[keep])
            report.windows_total += n_win
    finally:
        for reader in readers.values():
            reader.close()
    if not coords_start:
        raise ValueError("The window is larger than every selected chromosome")
    table = np.column_stack([np.concatenate(col) for col in columns])
    chrom_arr = np.concatenate(coords_chrom)
    start_arr = np.concatenate(coords_start)
    table, (chrom_arr, start_arr) = _filter(
        table, [chrom_arr, start_arr], skip_zeros, max_rows, seed, report
    )
    report.seconds = round(time.time() - began, 2)
    summary = RegionSummary(
        regions=tuple(f"{c}:{s}-{s + size}" for c, s in zip(chrom_arr, start_arr, strict=True)),
        samples=tuple(sources),
        chroms=tuple(chrom_arr.tolist()),
        starts=tuple(start_arr.tolist()),
        ends=tuple((start_arr + size).tolist()),
        values=table,
    )
    return summary, report


def _filter(table, arrays, skip_zeros, max_rows, seed, report):
    """Drop all-zero rows and optionally subsample; `arrays` are row-aligned annotations."""
    if skip_zeros and table.size:
        keep = np.any(table != 0, axis=1)
        report.zero_removed = int((~keep).sum())
        table, arrays = table[keep], [a[keep] for a in arrays]
    if max_rows and table.shape[0] > max_rows:
        pick = np.sort(np.random.default_rng(seed).choice(table.shape[0], max_rows, replace=False))
        table, arrays = table[pick], [a[pick] for a in arrays]
        report.subsampled = True
    report.windows_kept = int(table.shape[0])
    return table, arrays


def _alias(reader, chrom):
    if chrom in reader.sizes:
        return chrom
    key = chrom.lower().removeprefix("chr")
    return next((c for c in reader.sizes if c.lower().removeprefix("chr") == key), chrom)


def region_summary(
    sources: dict[str, str],
    bed_paths: list[str | Path],
    *,
    stat: str = "sum",
    skip_zeros: bool = False,
    max_rows: int | None = None,
    seed: int = 7201,
) -> tuple[RegionSummary, BinReport]:
    """Signal per BED region for every sample (multiBigwigSummary BED-file equivalent)."""
    began = time.time()
    if stat not in STATS:
        raise ValueError(f"Statistic must be one of {', '.join(STATS)}")
    features = [f for path in bed_paths for f in read_bed(path)]
    if not features:
        raise ValueError("The selected region files are empty")
    names = np.array([f.name or f"r{i}" for i, f in enumerate(features)], dtype=object)
    chroms = np.array([f.interval.chrom for f in features], dtype=object)
    starts = np.array([f.interval.start for f in features], dtype=np.int64)
    ends = np.array([f.interval.end for f in features], dtype=np.int64)
    table = np.zeros((len(features), len(sources)))
    report = BinReport(chroms=sorted(set(chroms.tolist())), windows_total=len(features))
    for column, path in enumerate(sources.values()):
        reader = SignalReader(path)
        try:
            if reader.dense:
                # Per-region reads: never materialise a whole chromosome of 1-bp values.
                widths = np.maximum(ends - starts, 1)
                total, cover, peak = (
                    np.zeros(len(starts)),
                    np.zeros(len(starts)),
                    np.full(len(starts), -np.inf),
                )
                for row, (chrom, a, b) in enumerate(zip(chroms, starts, ends, strict=True)):
                    block = reader.values(_alias(reader, chrom), a, b)
                    finite = np.isfinite(block)
                    total[row], cover[row] = block[finite].sum(), finite.sum()
                    if stat == "max" and finite.any():
                        peak[row] = block[finite].max()
                table[:, column] = _finish(stat, total, cover, peak, widths)
                continue
            for chrom in report.chroms:
                rows = np.flatnonzero(chroms == chrom)
                lo, hi = int(starts[rows].min()), int(ends[rows].max())
                s, e, v = reader.intervals(chrom, lo, hi)
                if not s.size and chrom not in reader.sizes:
                    alias = next(
                        (
                            c
                            for c in reader.sizes
                            if c.lower().removeprefix("chr") == chrom.lower().removeprefix("chr")
                        ),
                        None,
                    )
                    if alias:
                        s, e, v = reader.intervals(alias, lo, hi)
                f_sum_s, f_cov_s = _integrals(s, e, v, starts[rows].astype(float))
                f_sum_e, f_cov_e = _integrals(s, e, v, ends[rows].astype(float))
                total, cover = f_sum_e - f_sum_s, f_cov_e - f_cov_s
                peak = None
                if stat == "max":
                    peak = np.array(
                        [
                            v[(e > a) & (s < b)].max(initial=-np.inf)
                            for a, b in zip(starts[rows], ends[rows], strict=True)
                        ]
                    )
                table[rows, column] = _finish(
                    stat, total, cover, peak, np.maximum(ends[rows] - starts[rows], 1)
                )
        finally:
            reader.close()
    table, (names, chroms, starts, ends) = _filter(
        table, [names, chroms, starts, ends], skip_zeros, max_rows, seed, report
    )
    report.seconds = round(time.time() - began, 2)
    summary = RegionSummary(
        regions=tuple(names.tolist()),
        samples=tuple(sources),
        chroms=tuple(chroms.tolist()),
        starts=tuple(int(x) for x in starts),
        ends=tuple(int(x) for x in ends),
        values=table,
    )
    return summary, report


def window_matrix(path: str | Path, windows, n_bins: int, stat: str = "mean") -> np.ndarray:
    """Rows of binned signal for many windows (heatmap / profile matrices).

    Same semantics as per-window ``pyBigWig.stats(nBins=…)``: ``mean`` is over covered bases and
    0 where a bin has no data; bins beyond a chromosome edge are NaN. Windows may differ in
    length (scale-regions); bin edges are then fractional, exactly like pyBigWig.
    """
    if stat not in ("mean", "max", "sum"):
        raise ValueError("window_matrix supports mean, max or sum")
    reader = SignalReader(path)
    out = np.full((len(windows), n_bins), np.nan)
    try:
        by_chrom: dict[str, list[int]] = {}
        for index, (chrom, _start, _end) in enumerate(windows):
            by_chrom.setdefault(_alias(reader, chrom), []).append(index)
        frac = np.arange(n_bins + 1) / n_bins
        for chrom, rows in by_chrom.items():
            length = reader.sizes.get(chrom)
            if length is None:
                continue
            starts = np.array([windows[i][1] for i in rows], dtype=float)
            ends = np.array([windows[i][2] for i in rows], dtype=float)
            edges = starts[:, None] + (ends - starts)[:, None] * frac[None, :]
            inside = (edges[:, :-1] >= 0) & (edges[:, 1:] <= length)
            clipped = np.clip(edges, 0, length)
            if reader.dense or stat == "max":
                for k, row in enumerate(rows):
                    lo, hi = int(clipped[k, 0]), int(np.ceil(clipped[k, -1]))
                    if hi <= lo:
                        continue
                    block = reader.values(chrom, lo, hi)
                    finite = np.isfinite(block)
                    cut = np.clip(np.round(clipped[k] - lo).astype(np.int64), 0, hi - lo)
                    empty = cut[1:] <= cut[:-1]  # zero-width bins, e.g. beyond a chromosome end
                    starts_ = np.minimum(cut[:-1], hi - lo - 1)  # reduceat needs in-range indices
                    total = np.add.reduceat(np.where(finite, block, 0.0), starts_)
                    cover = np.add.reduceat(finite.astype(float), starts_)
                    total[empty], cover[empty] = 0.0, 0.0
                    if stat == "max":
                        peak = np.maximum.reduceat(np.where(finite, block, -np.inf), starts_)
                        peak[empty] = -np.inf
                        values = np.where(np.isfinite(peak), peak, 0.0)
                    elif stat == "sum":
                        values = total
                    else:
                        values = np.where(cover > 0, total / np.maximum(cover, 1e-12), 0.0)
                    out[row] = np.where(inside[k], values, np.nan)
                continue
            lo, hi = int(max(0, starts.min())), int(min(length, ends.max()))
            s, e, v = reader.intervals(chrom, lo, hi)
            f_sum, f_cov = _integrals(s, e, v, clipped.ravel())
            total = np.diff(f_sum.reshape(clipped.shape), axis=1)
            cover = np.diff(f_cov.reshape(clipped.shape), axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                values = (
                    total
                    if stat == "sum"
                    else np.where(cover > 0, total / np.maximum(cover, 1e-12), 0.0)
                )
            out[rows] = np.where(inside, values, np.nan)
    finally:
        reader.close()
    return out
