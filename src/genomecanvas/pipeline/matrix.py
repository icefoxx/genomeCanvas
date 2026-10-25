"""Grouped enrichment matrices for the heatmap pipeline."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from ..io.bed import read_bed
from ..io.bigwig import Engine, binned_bigwig, resolve_engine
from ..io.signal import infer_signal_format, query_signal
from ..model.feature import Feature
from ..profile.extract import extract_profile
from .metadata import SampleMeta

RefPoint = Literal["center", "TSS", "TES"]
SortUsing = Literal["mean", "median", "max", "min", "sum"]


@dataclass
class RegionCoord:
    """BED-like coordinates stored with an enrichment row."""

    chrom: str
    start: int
    end: int
    name: str
    score: str
    strand: str
    group: str


@dataclass
class EnrichmentMatrix:
    """Sample-concatenated region-by-bin matrix, analogous to computeMatrix."""

    values: np.ndarray  # type: ignore[type-arg]
    coords: list[RegionCoord]
    sample_labels: tuple[str, ...]
    sample_paths: tuple[str, ...]
    sample_boundaries: tuple[int, ...]
    group_labels: tuple[str, ...]
    group_boundaries: tuple[int, ...]
    upstream: int
    downstream: int
    body: int
    bin_size: int
    ref_point: str
    bins_per_sample: int
    meta_groups: tuple[str, ...]
    conditions: tuple[str, ...]
    colormaps: tuple[str, ...]
    profile_colors: tuple[str, ...]

    def sample_block(self, index: int) -> np.ndarray:  # type: ignore[type-arg]
        """Return the bins belonging to one sample."""
        start = self.sample_boundaries[index]
        end = self.sample_boundaries[index + 1]
        return self.values[:, start:end]

    def mean_profiles(self) -> dict[str, np.ndarray]:  # type: ignore[type-arg]
        """Return the column-mean profile for each sample."""
        return {
            label: np.nanmean(self.sample_block(index), axis=0)
            for index, label in enumerate(self.sample_labels)
        }

    def axis_positions(self) -> np.ndarray:  # type: ignore[type-arg]
        """Return bin-start coordinates in bp relative to the reference point."""
        if self.body:
            up = self.upstream
            body = self.body
            down = self.downstream
            left = np.linspace(-up, 0, max(1, up // self.bin_size), endpoint=False)
            mid = np.linspace(0, body, max(1, body // self.bin_size), endpoint=False)
            right = np.linspace(body, body + down, max(1, down // self.bin_size), endpoint=False)
            return np.concatenate([left, mid, right])[: self.bins_per_sample]
        half = (self.bins_per_sample * self.bin_size) // 2
        return np.arange(-half, half, self.bin_size, dtype=float)[: self.bins_per_sample]


def compute_region_matrix(
    samples: list[SampleMeta],
    beds: list[tuple[str, str]],
    *,
    flank: int = 2000,
    binsize: int = 50,
    ref_point: RefPoint = "center",
    body_length: int = 0,
    skip_zeros: bool = True,
    engine: Engine = "auto",
) -> EnrichmentMatrix:
    """Extract a grouped reference-point or scale-regions matrix."""
    if not samples:
        raise ValueError("at least one BigWig sample is required")
    if not beds:
        raise ValueError("at least one BED group is required")
    start_from: Literal["center", "start", "end"] = {
        "center": "center", "TSS": "start", "TES": "end"
    }[ref_point]
    sources = {item.label: item.path for item in samples}
    if body_length > 0 and ref_point in {"TSS", "TES"}:
        return _scale_regions_matrix(
            samples, beds, flank=flank, binsize=binsize, body_length=body_length,
            engine=engine, skip_zeros=skip_zeros,
        )
    group_labels: list[str] = []
    group_boundaries = [0]
    per_group: list[list[np.ndarray]] = []  # type: ignore[type-arg]
    coords: list[RegionCoord] = []
    for bed_path, group_label in beds:
        features = read_bed(bed_path)
        matrices = extract_profile(
            sources,
            bed=bed_path,
            start_from=start_from,
            up=flank,
            down=flank,
            binsize=binsize,
            engine=engine,
        )
        n_rows = matrices[0].values.shape[0] if matrices else 0
        if n_rows != len(features):
            n_rows = min(n_rows, len(features))
            features = features[:n_rows]
            matrices = [
                type(item)(
                    item.sample, item.values[:n_rows], item.up, item.down, item.binsize,
                    item.start_from, item.region_ids[:n_rows], item.source,
                )
                for item in matrices
            ]
        per_group.append([item.values for item in matrices])
        group_labels.append(group_label)
        group_boundaries.append(group_boundaries[-1] + n_rows)
        for feature in features:
            coords.append(
                RegionCoord(
                    feature.interval.chrom,
                    feature.interval.start,
                    feature.interval.end,
                    feature.name or ".",
                    str(feature.score if feature.score is not None else "."),
                    feature.interval.strand or ".",
                    group_label,
                )
            )
    stacked_per_sample = [
        np.vstack([block[index] for block in per_group]) if per_group else np.zeros((0, 1))
        for index in range(len(samples))
    ]
    values = np.hstack(stacked_per_sample) if stacked_per_sample else np.zeros((0, 0))
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    if skip_zeros and values.size:
        keep = np.any(values != 0, axis=1)
        values, coords, group_boundaries = _filter_rows(values, coords, group_boundaries, keep)
    bins_per_sample = int(stacked_per_sample[0].shape[1]) if stacked_per_sample else 0
    boundaries = [0]
    for _ in samples:
        boundaries.append(boundaries[-1] + bins_per_sample)
    return EnrichmentMatrix(
        values=np.asarray(values, dtype=np.float32),
        coords=coords,
        sample_labels=tuple(item.label for item in samples),
        sample_paths=tuple(item.path for item in samples),
        sample_boundaries=tuple(boundaries),
        group_labels=tuple(group_labels),
        group_boundaries=tuple(group_boundaries),
        upstream=flank,
        downstream=flank,
        body=0,
        bin_size=binsize,
        ref_point=ref_point,
        bins_per_sample=bins_per_sample,
        meta_groups=tuple(item.meta_group for item in samples),
        conditions=tuple(item.condition for item in samples),
        colormaps=tuple(item.colormap for item in samples),
        profile_colors=tuple(item.profile_color for item in samples),
    )


def _scale_regions_matrix(
    samples: list[SampleMeta],
    beds: list[tuple[str, str]],
    *,
    flank: int,
    binsize: int,
    body_length: int,
    engine: Engine,
    skip_zeros: bool,
) -> EnrichmentMatrix:
    n_up = max(1, flank // binsize) if flank else 0
    n_body = max(1, body_length // binsize)
    n_down = n_up
    bins_per_sample = n_up + n_body + n_down
    group_labels: list[str] = []
    group_boundaries = [0]
    coords: list[RegionCoord] = []
    row_blocks: list[np.ndarray] = []  # type: ignore[type-arg]
    resolved = resolve_engine(engine, batch=False)
    for bed_path, group_label in beds:
        features = read_bed(bed_path)
        group_rows = [
            _one_scale_region(item.path, features, n_up, n_body, n_down, flank, binsize, resolved)
            for item in samples
        ]
        stacked = np.hstack(group_rows) if group_rows else np.zeros((len(features), 0))
        row_blocks.append(stacked)
        group_labels.append(group_label)
        group_boundaries.append(group_boundaries[-1] + len(features))
        for feature in features:
            coords.append(
                RegionCoord(
                    feature.interval.chrom, feature.interval.start, feature.interval.end,
                    feature.name or ".", str(feature.score if feature.score is not None else "."),
                    feature.interval.strand or ".", group_label,
                )
            )
    values = np.vstack(row_blocks) if row_blocks else np.zeros((0, bins_per_sample * len(samples)))
    values = np.nan_to_num(values, nan=0.0)
    if skip_zeros and values.size:
        keep = np.any(values != 0, axis=1)
        values, coords, group_boundaries = _filter_rows(values, coords, group_boundaries, keep)
    boundaries = [0]
    for _ in samples:
        boundaries.append(boundaries[-1] + bins_per_sample)
    return EnrichmentMatrix(
        values=np.asarray(values, dtype=np.float32),
        coords=coords,
        sample_labels=tuple(item.label for item in samples),
        sample_paths=tuple(item.path for item in samples),
        sample_boundaries=tuple(boundaries),
        group_labels=tuple(group_labels),
        group_boundaries=tuple(group_boundaries),
        upstream=flank,
        downstream=flank,
        body=body_length,
        bin_size=binsize,
        ref_point="scale-regions",
        bins_per_sample=bins_per_sample,
        meta_groups=tuple(item.meta_group for item in samples),
        conditions=tuple(item.condition for item in samples),
        colormaps=tuple(item.colormap for item in samples),
        profile_colors=tuple(item.profile_color for item in samples),
    )


def _one_scale_region(
    path: str,
    features: list[Feature],
    n_up: int,
    n_body: int,
    n_down: int,
    flank: int,
    binsize: int,
    engine: str,
) -> np.ndarray:  # type: ignore[type-arg]
    kind = infer_signal_format(path)
    rows: list[np.ndarray] = []  # type: ignore[type-arg]
    for feature in features:
        interval = feature.interval
        plus = interval.strand != "-"
        tss = interval.start if plus else interval.end
        tes = interval.end if plus else interval.start
        if plus:
            five_start, five_end = max(0, tss - flank), tss
            three_start, three_end = tes, tes + flank
            body_start, body_end = interval.start, interval.end
            reverse_flanks = False
        else:
            five_start, five_end = tes, tes + flank
            three_start, three_end = max(0, tss - flank), tss
            body_start, body_end = interval.start, interval.end
            reverse_flanks = True
        parts = []
        for start, end, n_bins, reverse in (
            (five_start, max(five_start + 1, five_end), n_up, reverse_flanks),
            (body_start, max(body_start + 1, body_end), n_body, not plus),
            (three_start, max(three_start + 1, three_end), n_down, reverse_flanks),
        ):
            if kind == "bw":
                chunk = binned_bigwig(
                    path, interval.chrom, start, end, n_bins=n_bins, engine=engine  # type: ignore[arg-type]
                )
            else:
                chunk = query_signal(
                    path, interval.chrom, start, end, pixels=n_bins, engine=engine  # type: ignore[arg-type]
                ).values
            if reverse:
                chunk = chunk[::-1]
            parts.append(chunk)
        rows.append(np.concatenate(parts))
    if not rows:
        return np.zeros((0, n_up + n_body + n_down), dtype=float)
    return np.vstack(rows)


def _filter_rows(
    values: np.ndarray,  # type: ignore[type-arg]
    coords: list[RegionCoord],
    group_boundaries: list[int],
    keep: np.ndarray,  # type: ignore[type-arg]
) -> tuple[np.ndarray, list[RegionCoord], list[int]]:  # type: ignore[type-arg]
    kept_coords = [coord for coord, flag in zip(coords, keep, strict=True) if flag]
    new_bounds = [0]
    for left, right in zip(group_boundaries[:-1], group_boundaries[1:], strict=True):
        new_bounds.append(new_bounds[-1] + int(np.sum(keep[left:right])))
    return values[keep], kept_coords, new_bounds


def sort_matrix(
    matrix: EnrichmentMatrix,
    *,
    method: SortUsing = "mean",
    descending: bool = True,
    sample_indices: list[int] | None = None,
) -> EnrichmentMatrix:
    """Sort regions within each group using selected sample blocks."""
    if matrix.values.size == 0:
        return matrix
    indices = sample_indices or list(range(len(matrix.sample_labels)))
    blocks = [matrix.sample_block(index) for index in indices]
    joined = np.hstack(blocks) if blocks else matrix.values
    if method == "median":
        scores = np.nanmedian(joined, axis=1)
    elif method == "max":
        scores = np.nanmax(joined, axis=1)
    elif method == "min":
        scores = np.nanmin(joined, axis=1)
    elif method == "sum":
        scores = np.nansum(joined, axis=1)
    else:
        scores = np.nanmean(joined, axis=1)
    pieces: list[np.ndarray] = []  # type: ignore[type-arg]
    coords_out: list[RegionCoord] = []
    bounds = [0]
    for left, right in zip(matrix.group_boundaries[:-1], matrix.group_boundaries[1:], strict=True):
        group_idx = np.arange(left, right)
        ranked = group_idx[np.argsort(-scores[group_idx] if descending else scores[group_idx])]
        pieces.append(matrix.values[ranked])
        coords_out.extend(matrix.coords[int(index)] for index in ranked)
        bounds.append(bounds[-1] + len(ranked))
    values = np.vstack(pieces) if pieces else matrix.values
    return EnrichmentMatrix(
        values=values,
        coords=coords_out,
        sample_labels=matrix.sample_labels,
        sample_paths=matrix.sample_paths,
        sample_boundaries=matrix.sample_boundaries,
        group_labels=matrix.group_labels,
        group_boundaries=tuple(bounds),
        upstream=matrix.upstream,
        downstream=matrix.downstream,
        body=matrix.body,
        bin_size=matrix.bin_size,
        ref_point=matrix.ref_point,
        bins_per_sample=matrix.bins_per_sample,
        meta_groups=matrix.meta_groups,
        conditions=matrix.conditions,
        colormaps=matrix.colormaps,
        profile_colors=matrix.profile_colors,
    )


def save_matrix(matrix: EnrichmentMatrix, path: str | Path) -> Path:
    """Write values to ``.npz`` and metadata to a sibling ``.json``."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        destination,
        values=matrix.values,
        chrom=np.array([item.chrom for item in matrix.coords]),
        start=np.array([item.start for item in matrix.coords]),
        end=np.array([item.end for item in matrix.coords]),
        name=np.array([item.name for item in matrix.coords]),
        score=np.array([item.score for item in matrix.coords]),
        strand=np.array([item.strand for item in matrix.coords]),
        group=np.array([item.group for item in matrix.coords]),
    )
    header = {
        "sample_labels": list(matrix.sample_labels),
        "sample_paths": list(matrix.sample_paths),
        "sample_boundaries": list(matrix.sample_boundaries),
        "group_labels": list(matrix.group_labels),
        "group_boundaries": list(matrix.group_boundaries),
        "upstream": matrix.upstream,
        "downstream": matrix.downstream,
        "body": matrix.body,
        "bin_size": matrix.bin_size,
        "ref_point": matrix.ref_point,
        "bins_per_sample": matrix.bins_per_sample,
        "meta_groups": list(matrix.meta_groups),
        "conditions": list(matrix.conditions),
        "colormaps": list(matrix.colormaps),
        "profile_colors": list(matrix.profile_colors),
    }
    destination.with_suffix(".json").write_text(
        json.dumps(header, indent=2) + "\n", encoding="utf-8"
    )
    return destination


def load_matrix(path: str | Path) -> EnrichmentMatrix:
    """Load a matrix saved by :func:`save_matrix`."""
    source = Path(path)
    payload = np.load(source, allow_pickle=True)
    header = json.loads(source.with_suffix(".json").read_text(encoding="utf-8"))
    coords = [
        RegionCoord(
            str(chrom), int(start), int(end), str(name), str(score), str(strand), str(group)
        )
        for chrom, start, end, name, score, strand, group in zip(
            payload["chrom"], payload["start"], payload["end"], payload["name"],
            payload["score"], payload["strand"], payload["group"], strict=True,
        )
    ]
    return EnrichmentMatrix(
        values=np.asarray(payload["values"], dtype=np.float32),
        coords=coords,
        sample_labels=tuple(header["sample_labels"]),
        sample_paths=tuple(header["sample_paths"]),
        sample_boundaries=tuple(header["sample_boundaries"]),
        group_labels=tuple(header["group_labels"]),
        group_boundaries=tuple(header["group_boundaries"]),
        upstream=int(header["upstream"]),
        downstream=int(header["downstream"]),
        body=int(header["body"]),
        bin_size=int(header["bin_size"]),
        ref_point=str(header["ref_point"]),
        bins_per_sample=int(header["bins_per_sample"]),
        meta_groups=tuple(header["meta_groups"]),
        conditions=tuple(header["conditions"]),
        colormaps=tuple(header["colormaps"]),
        profile_colors=tuple(header["profile_colors"]),
    )


def write_mean_profiles(matrix: EnrichmentMatrix, path: str | Path) -> Path:
    """Write per-sample mean profiles as a TSV indexed by bp."""
    destination = Path(path)
    positions = matrix.axis_positions()
    profiles = matrix.mean_profiles()
    header = ["idx", *matrix.sample_labels]
    lines = ["\t".join(header)]
    for row, position in enumerate(positions):
        values = [f"{profiles[label][row]:.6g}" if row < profiles[label].size else "0"
                  for label in matrix.sample_labels]
        lines.append("\t".join([str(int(position)), *values]))
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination
