"""Genome-bin signal matrix for many BigWigs, built in parallel with per-sample checkpoints.

Every sample is binned over the same spans (taken from a reference file's chromosome sizes), so
columns line up even when files differ slightly in which contigs they carry. A killed run resumes:
finished samples are kept as ``<checkpoint>/<sha1>.npy`` and skipped next time.
"""

from __future__ import annotations

import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from ..io.bed import read_bed
from .genome_bins import (
    DEFAULT_SKIP,
    SignalReader,
    _alias,
    _finish,
    bin_chromosome,
    select_chromosomes,
)


def plan_spans(reference: str, size: int, chroms=None, skip: str = DEFAULT_SKIP) -> dict:
    reader = SignalReader(reference)
    try:
        sizes = dict(reader.sizes)
    finally:
        reader.close()
    spans = {c: (sizes[c] // size) * size for c in select_chromosomes(sizes, chroms, skip)}
    spans = {c: n for c, n in spans.items() if n}
    if not spans:
        raise ValueError("No chromosomes left after the include/skip lists")
    return spans


def bin_sample(path: str, spans: dict, size: int, stat: str = "mean") -> np.ndarray:
    reader = SignalReader(path)
    try:
        parts = []
        for chrom, length in spans.items():
            name = _alias(reader, chrom)
            if name not in reader.sizes:
                parts.append(np.zeros(length // size))
                continue
            total, cover, peak = bin_chromosome(reader, name, 0, length, size, stat)
            parts.append(_finish(stat, total, cover, peak, size))
        return np.concatenate(parts).astype(np.float32)
    finally:
        reader.close()


def _task(name, path, spans, size, stat, checkpoint):
    began = time.time()
    target = checkpoint / (
        hashlib.sha1(f"{path}|{size}|{stat}|{sorted(spans.items())}".encode()).hexdigest() + ".npy"
    )
    if target.is_file():
        return name, np.load(target), 0.0, True
    values = bin_sample(path, spans, size, stat)
    np.save(target, values)
    return name, values, time.time() - began, False


def build_matrix(
    sources: dict[str, str],
    out: str | Path,
    *,
    size: int = 10_000,
    stat: str = "mean",
    chroms=None,
    skip: str = DEFAULT_SKIP,
    blacklist: str | None = None,
    workers: int = 4,
    checkpoint: str | Path | None = None,
    log=print,
) -> dict:
    """Write ``out`` (.npz: values bins x samples float32, samples, chroms, starts, bin_size)."""
    out = Path(out)
    checkpoint = Path(checkpoint or out.with_suffix(".parts"))
    checkpoint.mkdir(parents=True, exist_ok=True)
    missing = {n: p for n, p in sources.items() if not Path(p).is_file()}
    usable = {n: p for n, p in sources.items() if n not in missing}
    if not usable:
        raise ValueError("None of the signal files exist")
    spans = plan_spans(next(iter(usable.values())), size, chroms, skip)
    columns, failed, began = {}, {}, time.time()
    with ProcessPoolExecutor(max_workers=max(1, workers)) as pool:
        jobs = {
            pool.submit(_task, n, p, spans, size, stat, checkpoint): n for n, p in usable.items()
        }
        for done, future in enumerate(as_completed(jobs), 1):
            name = jobs[future]
            try:
                _, values, seconds, cached = future.result()
                columns[name] = values
                log(f"[{done}/{len(jobs)}] {name}: {'checkpoint' if cached else f'{seconds:.1f}s'}")
            except Exception as exc:  # noqa: BLE001 - one unreadable file must not stop the run
                failed[name] = str(exc)
                log(f"[{done}/{len(jobs)}] {name}: FAILED {exc}")
    names = [n for n in usable if n in columns]
    if not names:
        first = next(iter(failed.items()), ("?", "unknown error"))
        raise ValueError(f"Every sample failed (e.g. {first[0]}: {first[1]})")
    values = np.column_stack([columns[n] for n in names])
    chrom_arr = np.concatenate([np.full(n // size, c) for c, n in spans.items()])
    start_arr = np.concatenate([np.arange(0, n, size, dtype=np.int64) for n in spans.values()])
    keep = np.ones(len(start_arr), dtype=bool)
    if blacklist:
        offset, positions = 0, {}
        for c, n in spans.items():
            positions[c] = offset
            offset += n // size
        for feature in read_bed(blacklist):
            c = feature.interval.chrom
            if c in positions:
                lo = feature.interval.start // size
                hi = min((feature.interval.end - 1) // size + 1, spans[c] // size)
                keep[positions[c] + lo : positions[c] + hi] = False
    values, chrom_arr, start_arr = values[keep], chrom_arr[keep], start_arr[keep]
    np.savez_compressed(
        out,
        values=values,
        samples=np.array(names),
        chroms=chrom_arr,
        starts=start_arr,
        bin_size=np.int64(size),
    )
    report = dict(
        out=str(out),
        samples=len(names),
        bins=int(values.shape[0]),
        blacklisted=int((~keep).sum()),
        missing_files=missing,
        failed=failed,
        seconds=round(time.time() - began, 1),
        bin_size=size,
        stat=stat,
    )
    out.with_suffix(".json").write_text(json.dumps(report, indent=2))
    return report
