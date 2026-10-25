"""BigWig adapters with selectable pyBigWig, pybigtools, and bwtool engines."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Literal

import math

import numpy as np

from ..core.coordinates import chrom_candidates

Aggregation = Literal["mean", "max", "min", "sum"]
Engine = Literal["auto", "pybigwig", "pybigtools", "bwtool"]
ResolvedEngine = Literal["pybigwig", "pybigtools", "bwtool"]

_ENGINE_OVERRIDE: Engine | None = None
_ENV_KEY = "GENOMECANVAS_BIGWIG_ENGINE"


def bwtool_available() -> bool:
    """Return whether the bwtool binary is on PATH."""
    return shutil.which("bwtool") is not None


def pybigwig_available() -> bool:
    """Return whether pyBigWig can be imported."""
    try:
        import pyBigWig  # type: ignore[import-not-found]  # noqa: F401
    except ImportError:
        return False
    return True


def pybigtools_available() -> bool:
    """Return whether pybigtools can be imported."""
    try:
        import pybigtools  # type: ignore[import-not-found]  # noqa: F401
    except ImportError:
        return False
    return True


def set_bigwig_engine(engine: Engine | None) -> None:
    """Set a process-wide engine override. Pass ``None`` to clear it."""
    global _ENGINE_OVERRIDE
    if engine is not None:
        _validate_engine_name(engine)
    _ENGINE_OVERRIDE = engine


def get_bigwig_engine() -> Engine:
    """Return the current requested engine, including override and environment."""
    if _ENGINE_OVERRIDE is not None:
        return _ENGINE_OVERRIDE
    value = os.environ.get(_ENV_KEY, "auto").strip().lower() or "auto"
    aliases = {
        "pybigwig": "pybigwig",
        "pybigwigs": "pybigwig",
        "pybigtools": "pybigtools",
        "bwtool": "bwtool",
        "auto": "auto",
    }
    if value not in aliases:
        raise ValueError(f"unknown {_ENV_KEY}={value!r}; use auto, pybigwig, pybigtools, or bwtool")
    return aliases[value]  # type: ignore[return-value]


@contextmanager
def using_bigwig_engine(engine: Engine) -> Iterator[ResolvedEngine]:
    """Temporarily force a BigWig engine for the current thread/process."""
    previous = _ENGINE_OVERRIDE
    set_bigwig_engine(engine)
    try:
        yield resolve_engine(engine)
    finally:
        set_bigwig_engine(previous)


def resolve_engine(requested: Engine = "auto", *, batch: bool = False) -> ResolvedEngine:
    """Resolve ``auto`` into a concrete backend.

    Batch jobs (profile matrices, many BED intervals) prefer bwtool.
    Single-locus queries prefer in-process pybigtools/pyBigWig so a subprocess
    is not started for every window.
    """
    _validate_engine_name(requested)
    choice = requested if requested != "auto" else get_bigwig_engine()
    if choice == "auto":
        if batch and bwtool_available():
            return "bwtool"
        if pybigtools_available():
            return "pybigtools"
        if pybigwig_available():
            return "pybigwig"
        if bwtool_available():
            return "bwtool"
        raise RuntimeError(
            "no BigWig engine available; install pyBigWig/pybigtools or put bwtool on PATH"
        )
    if choice == "bwtool":
        if not bwtool_available():
            raise RuntimeError("bwtool is not installed or not on PATH")
        return "bwtool"
    if choice == "pybigtools":
        if not pybigtools_available():
            raise RuntimeError("pybigtools is not installed; pip install pybigtools")
        return "pybigtools"
    if not pybigwig_available():
        raise RuntimeError("pyBigWig is not installed; pip install pyBigWig")
    return "pybigwig"


def bigwig_backend() -> str:
    """Return the resolved in-process Python backend name."""
    return resolve_engine("auto", batch=False)


def _validate_engine_name(engine: str) -> None:
    if engine not in {"auto", "pybigwig", "pybigtools", "bwtool"}:
        raise ValueError("engine must be auto, pybigwig, pybigtools, or bwtool")


def _match_chrom(available: list[str] | dict[str, Any], chrom: str) -> str:
    names = list(available) if not isinstance(available, dict) else list(available)
    lookup = {name: name for name in names}
    for candidate in chrom_candidates(chrom):
        if candidate in lookup:
            return candidate
    raise KeyError(f"chromosome {chrom!r} not found in BigWig chroms {names[:8]}")


def bigwig_chroms(path: str) -> dict[str, int]:
    """Return chromosome names and sizes from a BigWig header."""
    if pybigwig_available():
        import pyBigWig  # type: ignore[import-not-found, unused-ignore]

        handle = pyBigWig.open(path)
        try:
            chroms = handle.chroms() or {}
            return {str(key): int(value) for key, value in chroms.items()}
        finally:
            handle.close()
    if pybigtools_available():
        import pybigtools  # type: ignore[import-not-found, unused-ignore]

        handle = pybigtools.open(path)
        try:
            raw = handle.chroms()
            if isinstance(raw, dict):
                return {str(key): int(value) for key, value in raw.items()}
            return {str(name): int(size) for name, size in raw}
        finally:
            close = getattr(handle, "close", None)
            if callable(close):
                close()
    raise RuntimeError("reading BigWig chromosome sizes requires pyBigWig or pybigtools")


def read_bigwig(path: str, chrom: str, start: int, end: int) -> np.ndarray:  # type: ignore[type-arg]
    """Query per-base BigWig values for ``[start, end)``."""
    values = binned_bigwig(path, chrom, start, end, n_bins=max(1, end - start), aggregation="mean")
    return np.asarray(values, dtype=float)


def binned_bigwig(
    path: str,
    chrom: str,
    start: int,
    end: int,
    *,
    n_bins: int,
    aggregation: Aggregation = "mean",
    engine: Engine = "auto",
) -> np.ndarray:  # type: ignore[type-arg]
    """Return ``n_bins`` aggregated values over ``[start, end)``."""
    if end <= start or n_bins < 1:
        raise ValueError("BigWig query requires end > start and n_bins >= 1")
    resolved = resolve_engine(engine, batch=False)
    if resolved == "bwtool":
        return _bwtool_bins(path, chrom, start, end, n_bins, aggregation)
    if resolved == "pybigtools":
        return _pybigtools_bins(path, chrom, start, end, n_bins, aggregation)
    return _pybigwig_bins(path, chrom, start, end, n_bins, aggregation)


def _pybigwig_bins(
    path: str, chrom: str, start: int, end: int, n_bins: int, aggregation: Aggregation
) -> np.ndarray:  # type: ignore[type-arg]
    import pyBigWig  # type: ignore[import-not-found, unused-ignore]

    handle = pyBigWig.open(path)
    try:
        chroms = handle.chroms() or {}
        matched = _match_chrom(chroms, chrom)
        stats = handle.stats(matched, start, end, type=aggregation, nBins=n_bins)
    finally:
        handle.close()
    return np.array([0.0 if item is None else float(item) for item in stats], dtype=float)


def _pybigwig_matrix(
    path: str,
    windows: Sequence[tuple[str, int, int]],
    n_bins: int,
    aggregation: Aggregation,
) -> np.ndarray:  # type: ignore[type-arg]
    if aggregation in ("mean", "max", "sum") and windows:
        # Vectorised reader: one pass per chromosome instead of one stats() call per window.
        from ..pipeline.genome_bins import window_matrix

        return window_matrix(path, list(windows), n_bins, aggregation)
    import pyBigWig  # type: ignore[import-not-found, unused-ignore]

    handle = pyBigWig.open(path)
    try:
        chroms = handle.chroms() or {}
        rows: list[list[float]] = []
        for chrom, start, end in windows:
            matched = _match_chrom(chroms, chrom)
            rows.append(_window_stats(handle, chroms, matched, start, end, n_bins, aggregation))
    finally:
        handle.close()
    if not rows:
        return np.zeros((0, n_bins), dtype=float)
    return np.asarray(rows, dtype=float)


def _window_stats(handle, chroms, chrom, start, end, n_bins, aggregation) -> list[float]:
    """Binned stats for one window; bins outside the chromosome are NaN instead of an error."""
    length = chroms.get(chrom)
    if length is None:
        return [float("nan")] * n_bins
    if start >= 0 and end <= length:
        stats = handle.stats(chrom, start, end, type=aggregation, nBins=n_bins) or []
        return [0.0 if item is None else float(item) for item in stats]
    width = (end - start) / n_bins
    first = max(0, math.ceil((0 - start) / width))
    last = min(n_bins, math.floor((length - start) / width))
    row = [float("nan")] * n_bins
    if last > first:
        lo, hi = int(round(start + first * width)), int(round(start + last * width))
        stats = handle.stats(chrom, lo, min(hi, length), type=aggregation, nBins=last - first) or []
        row[first:last] = [0.0 if item is None else float(item) for item in stats]
    return row


def _pybigtools_bins(
    path: str, chrom: str, start: int, end: int, n_bins: int, aggregation: Aggregation
) -> np.ndarray:  # type: ignore[type-arg]
    import pybigtools  # type: ignore[import-not-found, unused-ignore]

    handle = pybigtools.open(path)
    try:
        chroms = list(handle.chroms())
        matched = _match_chrom(chroms, chrom)
        raw = np.asarray(handle.values(matched, start, end), dtype=float)
    finally:
        close = getattr(handle, "close", None)
        if callable(close):
            close()
    raw = np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0)
    if raw.size == 0:
        return np.zeros(n_bins, dtype=float)
    if raw.size == n_bins:
        return raw
    edges = np.linspace(0, raw.size, n_bins + 1, dtype=int)
    out = np.zeros(n_bins, dtype=float)
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        chunk = raw[left:right]
        if not chunk.size:
            continue
        if aggregation == "max":
            out[index] = float(np.max(chunk))
        elif aggregation == "min":
            out[index] = float(np.min(chunk))
        elif aggregation == "sum":
            out[index] = float(np.sum(chunk))
        else:
            out[index] = float(np.mean(chunk))
    return out


def _bwtool_float(value: str | None) -> float:
    """Parse a bwtool numeric field; uncovered intervals are written as NA."""
    if value is None:
        return 0.0
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", "NULL", "."}:
        return 0.0
    try:
        number = float(text)
    except ValueError:
        return 0.0
    if not np.isfinite(number):
        return 0.0
    return number


def _rewrite_bed_chroms_for_bigwig(bed: Path, bigwig: str, destination: Path) -> Path:
    """Rewrite BED chromosome names to match a BigWig header."""
    try:
        available = list(bigwig_chroms(bigwig))
    except RuntimeError:
        destination.write_bytes(Path(bed).read_bytes())
        return destination
    cache: dict[str, str] = {}
    lines: list[str] = []
    for raw in Path(bed).read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("track"):
            continue
        fields = raw.split("\t") if "\t" in raw else raw.split()
        if len(fields) < 3:
            continue
        chrom = fields[0]
        if chrom not in cache:
            try:
                cache[chrom] = _match_chrom(available, chrom)
            except KeyError:
                cache[chrom] = chrom
        fields[0] = cache[chrom]
        try:
            fields[1] = str(int(float(fields[1])))
            fields[2] = str(int(float(fields[2])))
        except (TypeError, ValueError):
            continue
        if len(fields) > 4:
            try:
                fields[4] = str(int(round(float(fields[4]))))
            except (TypeError, ValueError):
                fields[4] = "0"
        if len(fields) > 5 and fields[5] not in {"+", "-"}:
            fields[5] = "."
        lines.append("\t".join(fields[:6] if len(fields) > 6 else fields))
    destination.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return destination


def _run_bwtool(arguments: list[str]) -> None:
    try:
        subprocess.run(["bwtool", *arguments], check=True, capture_output=True, text=True)
    except FileNotFoundError as error:
        raise RuntimeError("bwtool is not installed or not on PATH") from error
    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or "").strip()
        raise RuntimeError(f"bwtool failed: {detail or error}") from error


def _bwtool_bins(
    path: str, chrom: str, start: int, end: int, n_bins: int, aggregation: Aggregation
) -> np.ndarray:  # type: ignore[type-arg]
    matched = chrom
    try:
        matched = _match_chrom(list(bigwig_chroms(path)), chrom)
    except (KeyError, RuntimeError):
        matched = chrom
    edges = np.linspace(start, end, n_bins + 1, dtype=int)
    with tempfile.TemporaryDirectory(prefix="gcanvas-bwtool-") as tmp:
        bed = Path(tmp) / "bins.bed"
        out = Path(tmp) / "summary.txt"
        lines = [
            f"{matched}\t{int(left)}\t{int(right)}\tbin{index}"
            for index, (left, right) in enumerate(zip(edges[:-1], edges[1:], strict=True))
            if right > left
        ]
        bed.write_text("\n".join(lines) + "\n", encoding="utf-8")
        _run_bwtool(
            [
                "summary",
                str(bed),
                path,
                str(out),
                "-keep-bed",
                "-header",
                "-skip-median",
                "-with-sum",
            ]
        )
        table = _read_bwtool_summary(out)
    column = {"mean": "mean", "max": "max", "min": "min", "sum": "sum"}[aggregation]
    values = np.array([_bwtool_float(row.get(column)) for row in table], dtype=float)
    if values.size < n_bins:
        values = np.pad(values, (0, n_bins - values.size))
    return values[:n_bins]


def _read_bwtool_summary(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    header = [item.lstrip("#") for item in lines[0].split()]
    rows: list[dict[str, str]] = []
    for raw in lines[1:]:
        fields = raw.split()
        if len(fields) < len(header):
            continue
        rows.append(dict(zip(header, fields, strict=False)))
    return rows


def bwtool_summary(path: str, chrom: str, start: int, end: int, *, binsize: int) -> np.ndarray:  # type: ignore[type-arg]
    """Extract binned means with bwtool when the binary is installed."""
    n_bins = max(1, (end - start + binsize - 1) // binsize)
    return binned_bigwig(
        path, chrom, start, end, n_bins=n_bins, aggregation="mean", engine="bwtool"
    )


def bwtool_matrix(
    bigwig: str,
    bed: Path,
    *,
    up: int,
    down: int,
    output: Path | None = None,
    binsize: int | None = None,
    start_from: Literal["start", "center", "end"] = "center",
) -> Path:
    """Run ``bwtool matrix`` and return the matrix path."""
    resolve_engine("bwtool")
    destination = output or Path(tempfile.mkstemp(suffix=".matrix")[1])
    with tempfile.TemporaryDirectory(prefix="gcanvas-bwtool-bed-") as tmp:
        remapped = Path(tmp) / "regions.bed"
        _rewrite_bed_chroms_for_bigwig(Path(bed), bigwig, remapped)
        arguments = ["matrix", f"{up}:{down}", str(remapped), bigwig, str(destination)]
        if start_from == "start":
            arguments.append("-starts")
        elif start_from == "end":
            arguments.append("-ends")
        if binsize is not None and binsize > 1:
            arguments.append(f"-tiled-averages={int(binsize)}")
        _run_bwtool(arguments)
    return Path(destination)


def read_bwtool_matrix(
    bigwig: str,
    bed: str | Path,
    *,
    up: int,
    down: int,
    binsize: int = 50,
    start_from: Literal["start", "center", "end"] = "center",
) -> np.ndarray:  # type: ignore[type-arg]
    """Return a region-by-bin matrix from ``bwtool matrix``."""
    n_bins = max(1, (up + down) // max(1, binsize))
    with tempfile.TemporaryDirectory(prefix="gcanvas-bwtool-") as tmp:
        out = Path(tmp) / "matrix.txt"
        bwtool_matrix(
            bigwig,
            Path(bed),
            up=up,
            down=down,
            output=out,
            binsize=binsize,
            start_from=start_from,
        )
        rows: list[list[float]] = []
        for raw in out.read_text(encoding="utf-8").splitlines():
            if not raw.strip():
                continue
            values = [_bwtool_float(item) for item in raw.split()]
            rows.append(values)
    if not rows:
        return np.zeros((0, n_bins), dtype=float)
    matrix = np.asarray(rows, dtype=float)
    if matrix.shape[1] < n_bins:
        matrix = np.pad(matrix, ((0, 0), (0, n_bins - matrix.shape[1])))
    return matrix[:, :n_bins]


def read_bwtool_region_sums(bigwig: str, bed: str | Path) -> np.ndarray:  # type: ignore[type-arg]
    """Return per-interval sums from ``bwtool summary -with-sum``."""
    resolve_engine("bwtool")
    with tempfile.TemporaryDirectory(prefix="gcanvas-bwtool-") as tmp:
        remapped = Path(tmp) / "regions.bed"
        _rewrite_bed_chroms_for_bigwig(Path(bed), bigwig, remapped)
        out = Path(tmp) / "summary.txt"
        _run_bwtool(
            [
                "summary",
                str(remapped),
                bigwig,
                str(out),
                "-keep-bed",
                "-header",
                "-skip-median",
                "-with-sum",
            ]
        )
        table = _read_bwtool_summary(out)
    return np.array([_bwtool_float(row.get("sum")) for row in table], dtype=float)


def python_profile_matrix(
    path: str,
    windows: Sequence[tuple[str, int, int]],
    *,
    n_bins: int,
    aggregation: Aggregation = "mean",
    engine: Engine = "auto",
) -> np.ndarray:  # type: ignore[type-arg]
    """Extract a profile matrix with an in-process BigWig library, one file open."""
    resolved = resolve_engine(engine, batch=False)
    if resolved == "pybigwig":
        return _pybigwig_matrix(path, windows, n_bins, aggregation)
    rows = [
        binned_bigwig(
            path, chrom, start, end, n_bins=n_bins, aggregation=aggregation, engine=resolved
        )
        for chrom, start, end in windows
    ]
    return np.vstack(rows) if rows else np.zeros((0, n_bins), dtype=float)
