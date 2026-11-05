"""Sample QC for large ChIP-seq collections: are each sample's mark and group labels plausible?

Input is a signal matrix (genome bins x samples, from ``gcanvas-sampleqc matrix`` or deepTools
``multiBigwigSummary --outRawCounts``) plus a sample sheet naming each sample's mark, group and
batch. Two layers:

1. Mark check over all samples: the k nearest neighbours in correlation space vote for a mark.
   A sample is flagged when the vote disagrees with its label.
2. Group check within each mark: a sample's correlation with the nearest other member of its own
   group is compared with its correlation with the nearest member of every other group. With 1-3
   samples per group, per-sample comparisons are the only ones that stay fair: a group centroid is
   smoother than any single sample (so it wins by construction), and a median over two members is
   pulled down by a single swapped sample, which would flag its innocent group-mates too.
   Samples judged low quality are left out of everyone else's comparisons.

Every sample gets one status, most severe first: ``mark_mismatch``, ``low_quality``, ``swap``,
``mislabel``, ``mislabel_batch`` (the nearest neighbour shares the batch, so this may be a batch
effect), ``undetermined`` (no other member of its group in this mark), ``ok``.
"""

from __future__ import annotations

import csv
import io
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

STATUS_ORDER = (
    "mark_mismatch",
    "low_quality",
    "swap",
    "mislabel",
    "mislabel_batch",
    "undetermined",
    "ok",
)
SIGNAL_EXT = re.compile(r"\.(bw|bigwig|bedgraph|bdg|bg|bam)$", re.IGNORECASE)
GUESS = {
    "sample": ("sample_id", "sample", "sample_name", "name", "id", "label"),
    "mark": ("mark", "target", "antibody_target", "histone", "factor", "chip", "assay"),
    "group": ("group", "condition", "model", "cell_line", "cell", "genotype", "tissue"),
    "batch": ("batch", "run", "run_id", "lane", "flowcell"),
}
MISSING = {"", "na", "n/a", "nan", "none", "null", "-"}


# ---------------------------------------------------------------- inputs


@dataclass
class SignalMatrix:
    samples: list[str]
    values: np.ndarray  # features x samples, float32
    chroms: np.ndarray | None = None
    starts: np.ndarray | None = None
    bin_size: int | None = None


def _clean_name(name: str) -> str:
    return str(name).strip().strip("'\"").lstrip("#").strip("'\"")


def matrix_samples(path: str | Path) -> list[str]:
    """Sample names of a matrix without loading its values."""
    path = Path(path)
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as data:
            return [str(s) for s in data["samples"]]
    with _open_text(path) as handle:
        header = handle.readline().rstrip("\n\r")
    names = [_clean_name(c) for c in re.split(r"\t|,", header)]
    return names[_coord_columns(names) :]


def _coord_columns(names: list[str]) -> int:
    lowered = [n.lower() for n in names[:3]]
    if len(lowered) == 3 and lowered[0] in ("chr", "chrom", "chromosome", "seqnames"):
        return 3
    return 0


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        import gzip

        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def load_matrix(path: str | Path) -> SignalMatrix:
    """Read a ``.npz`` written by ``gcanvas-sampleqc matrix`` or a tab/comma separated table
    whose first three columns are chrom, start, end (deepTools ``--outRawCounts`` included)."""
    path = Path(path)
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as data:
            values = np.asarray(data["values"], dtype=np.float32)
            samples = [str(s) for s in data["samples"]]
            chroms = data["chroms"] if "chroms" in data else None
            starts = data["starts"] if "starts" in data else None
            size = int(data["bin_size"]) if "bin_size" in data else None
        if values.shape[1] != len(samples):
            raise ValueError(f"{path.name}: {values.shape[1]} columns but {len(samples)} names")
        return SignalMatrix(samples, values, chroms, starts, size)
    with _open_text(path) as handle:
        header = handle.readline().rstrip("\n\r")
        delimiter = "\t" if "\t" in header else ","
        names = [_clean_name(c) for c in header.split(delimiter)]
        skip = _coord_columns(names)
        rows, chroms, starts = [], [], []
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n\r").split(delimiter)
            if skip:
                chroms.append(parts[0])
                starts.append(int(float(parts[1])))
            rows.append([float(v) if v not in ("", "nan", "NA") else np.nan for v in parts[skip:]])
    values = np.asarray(rows, dtype=np.float32)
    if values.ndim != 2 or values.shape[1] != len(names) - skip:
        raise ValueError(f"{path.name}: rows do not match the {len(names) - skip} sample columns")
    return SignalMatrix(
        names[skip:],
        np.nan_to_num(values, nan=0.0),
        np.asarray(chroms) if skip else None,
        np.asarray(starts) if skip else None,
    )


def read_sheet(path: str | Path) -> tuple[list[str], list[dict]]:
    """Columns and rows of a tab or comma separated sample sheet."""
    with _open_text(Path(path)) as handle:
        text = handle.read()
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("The sample sheet is empty")
    delimiter = "\t" if "\t" in lines[0] else ","
    reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter=delimiter)
    columns = [c.strip() for c in reader.fieldnames or []]
    rows = [{c.strip(): (v or "").strip() for c, v in row.items() if c} for row in reader]
    return columns, rows


def parse_filters(text: str) -> list[tuple[str, bool, set[str]]]:
    """``status=PASS; target!=INPUT,IgG`` -> [(column, keep_if_match, values)]."""
    filters = []
    for part in re.split(r"[;\n]", text or ""):
        part = part.strip()
        if not part:
            continue
        match = re.fullmatch(r"([^=!]+?)\s*(!=|=)\s*(.*)", part)
        if not match:
            raise ValueError(f"Filter {part!r} should look like column=value or column!=value")
        column, op, raw = match.groups()
        values = {v.strip() for v in raw.split(",")}
        filters.append((column.strip(), op == "=", values))
    return filters


def apply_filters(rows: list[dict], filters) -> list[dict]:
    out = rows
    for column, keep, values in filters:
        if out and column not in out[0]:
            raise ValueError(f"Filter column {column!r} is not in the sample sheet")
        out = [r for r in out if (r.get(column, "") in values) == keep]
    return out


def _norm(name: str) -> str:
    return SIGNAL_EXT.sub("", Path(_clean_name(name)).name)


def match_samples(rows: list[dict], column: str, samples: list[str]) -> dict[int, int]:
    """Sheet row index -> matrix column index, matching names with paths/extensions removed."""
    index = {}
    for j, name in enumerate(samples):
        index.setdefault(_norm(name), j)
    return {
        i: index[_norm(r.get(column, ""))]
        for i, r in enumerate(rows)
        if _norm(r.get(column, "")) in index
    }


def guess_columns(columns: list[str], rows: list[dict], samples: list[str] | None = None) -> dict:
    """Best-guess sample / mark / group / batch columns from names, cardinality and matches."""
    lower = {c.lower(): c for c in columns}
    n = max(1, len(rows))
    distinct = {c: len({r.get(c, "") for r in rows}) for c in columns}
    out: dict[str, str | None] = {}
    if samples:
        hits = {c: len(match_samples(rows, c, samples)) for c in columns}
        best = max(columns, key=lambda c: hits[c], default=None)
        out["sample"] = best if best and hits[best] else None
    for role, names in GUESS.items():
        if out.get(role):
            continue
        for name in names:
            column = lower.get(name)
            repeated = column and (role == "sample" or 1 < distinct[column] < n)
            if repeated and column not in out.values():
                out[role] = column
                break
        else:
            if role == "sample":
                out[role] = next((c for c in columns if distinct[c] == n), None)
            elif role == "batch":
                out[role] = next((c for c in columns if "batch" in c.lower()), None)
            else:
                out[role] = None
    return out


def preview(sheet: str | Path, matrix: str | Path | None = None, options: dict | None = None):
    """Everything the column-mapping dialog shows: rows, column stats, guesses, matches, tiers."""
    options = options or {}
    columns, rows = read_sheet(sheet)
    samples = matrix_samples(matrix) if matrix else None
    guesses = guess_columns(columns, rows, samples)
    stats = []
    for c in columns:
        counts = Counter(r.get(c, "") for r in rows)
        stats.append(dict(column=c, distinct=len(counts), top=counts.most_common(6)))
    chosen = {
        k: options.get(f"qc_{k}") or guesses.get(k) for k in ("sample", "mark", "group", "batch")
    }
    out = dict(
        columns=columns,
        rows=rows[:200],
        n_rows=len(rows),
        stats=stats,
        guesses=guesses,
        chosen=chosen,
        matrix_samples=len(samples) if samples is not None else None,
    )
    try:
        kept = apply_filters(rows, parse_filters(options.get("qc_filters", "")))
    except ValueError as exc:
        out["filter_error"] = str(exc)
        kept = rows
    out["n_kept"] = len(kept)
    if samples is not None and chosen["sample"]:
        matched = match_samples(kept, chosen["sample"], samples)
        out["n_matched"] = len(matched)
        out["unmatched_rows"] = [
            kept[i].get(chosen["sample"], "") for i in range(len(kept)) if i not in matched
        ][:20]
        kept = [kept[i] for i in matched]
    if chosen["mark"] and chosen["group"]:
        marks = Counter(r.get(chosen["mark"], "") for r in kept)
        table = []
        for mark, count in sorted(marks.items(), key=lambda kv: -kv[1]):
            groups = Counter(
                r.get(chosen["group"], "") for r in kept if r.get(chosen["mark"]) == mark
            )
            table.append(
                dict(
                    mark=mark,
                    samples=count,
                    groups=len(groups),
                    singletons=sum(1 for v in groups.values() if v == 1),
                    tier=tier(count),
                )
            )
        out["marks"] = table
    return out


def tier(n: int) -> str:
    return "S" if n <= 30 else "M" if n <= 100 else "L" if n <= 500 else "XL"


# ---------------------------------------------------------------- analysis


@dataclass
class QCResult:
    samples: list[str]
    marks: list[str]
    groups: list[str]
    batches: list[str]
    rows: list[dict]  # one per sample, flags table
    confusion: dict  # labelled mark -> predicted mark -> count
    mark_summary: list[dict]
    group_matrices: dict  # mark -> (group names, G x G medians, G x G pair counts)
    overview_pca: np.ndarray  # samples x 2
    overview_variance: np.ndarray
    per_mark: dict = field(default_factory=dict)  # mark -> dict(index, r, pca, variance)
    params: dict = field(default_factory=dict)


def _features(values: np.ndarray, top: int, transform: str) -> np.ndarray:
    x = np.log1p(np.clip(values, 0, None)) if transform == "log1p" else values.astype(np.float32)
    x = x[np.any(x != 0, axis=1)]
    if top and x.shape[0] > top:
        variance = x.var(axis=1)
        x = x[np.sort(np.argpartition(variance, -top)[-top:])]
    return np.asarray(x, dtype=np.float64)


def _correlate(x: np.ndarray, method: str) -> np.ndarray:
    if method == "spearman":
        from scipy.stats import rankdata

        x = rankdata(x, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.corrcoef(x, rowvar=False)
    return np.nan_to_num(np.atleast_2d(r), nan=0.0)


def _pca(x: np.ndarray, components: int = 2):
    table = x.T - x.T.mean(axis=0, keepdims=True)
    if table.shape[0] < 2:
        return np.zeros((table.shape[0], components)), np.zeros(components)
    _u, singular, vt = np.linalg.svd(table, full_matrices=False)
    coords = table @ vt[:components].T
    power = singular**2
    variance = power[:components] / power.sum() * 100 if power.sum() else np.zeros(components)
    if coords.shape[1] < components:
        pad = components - coords.shape[1]
        coords = np.pad(coords, ((0, 0), (0, pad)))
        variance = np.pad(variance, (0, pad))
    return coords, variance


def run_qc(
    matrix: SignalMatrix,
    meta: list[dict],
    *,
    k: int = 5,
    top: int = 5000,
    transform: str = "log1p",
    method: str = "pearson",
    z: float = 3.0,
    margin: float = 0.02,
    skip_group_marks: tuple[str, ...] = ("INPUT", "IgG"),
) -> QCResult:
    """`meta` rows (one per matrix column, same order) carry sample, mark, group, batch."""
    n = len(meta)
    if n < 3:
        raise ValueError("Sample QC needs at least three matched samples")
    samples = [m["sample"] for m in meta]
    marks = [m["mark"] or "(none)" for m in meta]
    groups = [m["group"] or "(none)" for m in meta]
    batches = [m.get("batch") or "" for m in meta]
    rows = [
        dict(
            sample=samples[i],
            mark=marks[i],
            group=groups[i],
            batch=batches[i],
            status="ok",
            reasons=[],
        )
        for i in range(n)
    ]

    # Layer 1: mark labels over all samples.
    x_all = _features(matrix.values, top, transform)
    r_all = _correlate(x_all, method)
    np.fill_diagonal(r_all, -np.inf)
    mark_counts = Counter(marks)
    confusion: dict = {}
    for i in range(n):
        k_eff = min(k, n - 1)
        nearest = np.argsort(-r_all[i])[:k_eff]
        votes = Counter(marks[j] for j in nearest)
        best = max(votes.values())
        tied = {m for m, c in votes.items() if c == best}
        predicted = next(marks[j] for j in nearest if marks[j] in tied)
        own = votes.get(marks[i], 0) / k_eff
        rows[i].update(predicted_mark=predicted, mark_vote=round(own, 2))
        confusion.setdefault(marks[i], Counter())[predicted] += 1
        if mark_counts[marks[i]] > 1 and predicted != marks[i] and own < 0.5:
            rows[i]["reasons"].append(
                f"mark_mismatch: neighbours vote {predicted} ({votes[predicted]}/{k_eff})"
            )
    overview, overview_var = _pca(x_all)

    # Layer 2: group labels within each mark.
    summary = []
    matrices = {}
    per_mark = {}
    skip = {s.lower() for s in skip_group_marks}
    for mark in sorted(mark_counts, key=lambda m: -mark_counts[m]):
        idx = np.array([i for i in range(n) if marks[i] == mark])
        info = dict(mark=mark, samples=len(idx), tier=tier(len(idx)))
        info["groups"] = len({groups[i] for i in idx})
        info["singletons"] = sum(1 for c in Counter(groups[i] for i in idx).values() if c == 1)
        if len(idx) < 3 or mark.lower() in skip:
            info["note"] = "skipped (control)" if mark.lower() in skip else "fewer than 3 samples"
            summary.append(info)
            continue
        x = _features(matrix.values[:, idx], top, transform)
        r = _correlate(x, method)
        off = r.copy()
        np.fill_diagonal(off, np.nan)
        nn = np.nanargmax(off, axis=1)
        max_r = off[np.arange(len(idx)), nn]
        local_groups = np.array([groups[i] for i in idx])
        names = sorted(set(local_groups.tolist()))
        members = {g: np.flatnonzero(local_groups == g) for g in names}
        cross = np.where(local_groups[:, None] != local_groups[None, :], off, np.nan)
        with np.errstate(all="ignore"):
            best_cross = np.nanmax(np.where(np.isnan(cross), -np.inf, cross), axis=1)
        # Low quality = unlike everything: low best r overall AND low best r to other groups.
        # The second test spares a sample from a distinctive group whose replicates agree, and
        # judges a group's only member by the same yardstick as everyone else.
        floor_any = _floor(max_r, z)
        floor_cross = (
            _floor(best_cross[np.isfinite(best_cross)], z) if np.isfinite(best_cross).any() else -1
        )
        poor = (max_r < floor_any) & ~(best_cross >= floor_cross)
        for a, i in enumerate(idx):
            row = rows[i]
            j = idx[nn[a]]
            row.update(
                max_r=round(float(max_r[a]), 4),
                nn_sample=samples[j],
                nn_group=groups[j],
                nn_batch=batches[j],
            )
            if poor[a]:
                row["reasons"].append(
                    f"low_quality: best r {max_r[a]:.3f} < {floor_any:.3f} and best r to other "
                    f"groups {best_cross[a]:.3f} < {floor_cross:.3f}"
                )
                continue
            # Nearest usable member per group; low-quality samples are left out so one bad
            # replicate cannot make its group-mates look mislabelled.
            nearest = {}
            for g in names:
                usable = members[g][(members[g] != a) & ~poor[members[g]]]
                if usable.size:
                    nearest[g] = float(r[a, usable].max())
            own = nearest.pop(groups[i], None)
            if nearest:
                best_group = max(nearest, key=nearest.get)
                row.update(best_other_group=best_group, best_other_r=round(nearest[best_group], 4))
            if own is None:
                alone = len(members[groups[i]]) == 1
                why = "only member" if alone else "other members are low quality"
                row["reasons"].append(f"undetermined: {why} of {groups[i]} in {mark}")
                continue
            row["own_r"] = round(own, 4)
            if nearest:
                row["margin"] = round(own - nearest[best_group], 4)
        for i in idx:
            row = rows[i]
            if row.get("margin") is None or row["margin"] >= -margin:
                continue
            partner = next(
                (
                    samples[idx[b]]
                    for b in range(len(idx))
                    if groups[idx[b]] == row["best_other_group"]
                    and rows[idx[b]].get("best_other_group") == groups[i]
                    and (rows[idx[b]].get("margin") or 0) < -margin
                ),
                None,
            )
            text = (
                f"{row['own_r']:.3f} to nearest {groups[i]} vs "
                f"{row['best_other_r']:.3f} to nearest {row['best_other_group']}"
            )
            if partner:
                row["swap_with"] = partner
                row["reasons"].append(f"swap: {text}; {partner} points back")
            elif batches[i] and row["nn_batch"] == batches[i] and row["nn_group"] != groups[i]:
                row["reasons"].append(f"mislabel_batch: {text}; nearest neighbour shares batch")
            else:
                row["reasons"].append(f"mislabel: {text}")
        within, between, same_batch, other_batch = [], [], [], []
        size = len(names)
        grid = np.full((size, size), np.nan)
        pairs = np.zeros((size, size), dtype=int)
        for p, g in enumerate(names):
            for q, h in enumerate(names):
                block = r[np.ix_(members[g], members[h])]
                values = block[np.triu_indices(len(members[g]), 1)] if p == q else block.ravel()
                pairs[p, q] = values.size
                if values.size:
                    grid[p, q] = np.median(values)
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                if local_groups[a] == local_groups[b]:
                    within.append(r[a, b])
                    continue
                between.append(r[a, b])
                if batches[idx[a]]:
                    target = same_batch if batches[idx[a]] == batches[idx[b]] else other_batch
                    target.append(r[a, b])
        info.update(
            median_r_within_group=_med(within),
            median_r_between_groups=_med(between),
            batch_delta=round(_med(same_batch) - _med(other_batch), 4)
            if same_batch and other_batch
            else None,
            low_quality_floor=round(float(floor_any), 4),
        )
        matrices[mark] = (names, grid, pairs)
        coords, variance = _pca(x)
        per_mark[mark] = dict(index=idx, r=r, pca=coords, variance=variance)
        summary.append(info)

    for row in rows:
        kinds = [reason.split(":", 1)[0] for reason in row["reasons"]]
        row["status"] = next((s for s in STATUS_ORDER if s in kinds), "ok")
        row["n_flags"] = sum(1 for kind in kinds if kind != "undetermined")
    for info in summary:
        info["flagged"] = sum(
            1
            for row in rows
            if row["mark"] == info["mark"] and row["status"] not in ("ok", "undetermined")
        )
    return QCResult(
        samples=samples,
        marks=marks,
        groups=groups,
        batches=batches,
        rows=rows,
        confusion={m: dict(c) for m, c in confusion.items()},
        mark_summary=summary,
        group_matrices=matrices,
        overview_pca=overview,
        overview_variance=overview_var,
        per_mark=per_mark,
        params=dict(k=k, top=top, transform=transform, method=method, z=z, margin=margin),
    )


def _floor(values, z: float) -> float:
    """Median - z robust SDs in Fisher-z units, back on the r scale. The spread has a floor of
    0.15 (z units) so near-identical replicates are not flagged over third-decimal differences."""
    zs = np.arctanh(np.clip(np.asarray(values, dtype=float), -0.999999, 0.999999))
    centre = float(np.median(zs))
    spread = max(1.4826 * float(np.median(np.abs(zs - centre))), 0.15)
    return float(np.tanh(centre - z * spread))


def _med(values) -> float | None:
    return round(float(np.median(values)), 4) if len(values) else None


FLAG_COLUMNS = (
    "sample", "mark", "group", "batch", "status", "n_flags", "predicted_mark", "mark_vote",
    "max_r", "own_r", "best_other_group", "best_other_r", "margin", "nn_sample", "nn_group",
    "nn_batch", "swap_with", "reasons",
)  # fmt: skip


def sorted_rows(result: QCResult) -> list[dict]:
    rank = {s: i for i, s in enumerate(STATUS_ORDER)}
    return sorted(result.rows, key=lambda r: (rank[r["status"]], r.get("margin") or 0, r["mark"]))


def write_tables(result: QCResult, destination: Path) -> None:
    def fmt(value):
        if value is None:
            return ""
        if isinstance(value, list):
            return " | ".join(value)
        return str(value)

    with open(destination / "qc_flags.tsv", "w", encoding="utf-8") as handle:
        handle.write("\t".join(FLAG_COLUMNS) + "\n")
        for row in sorted_rows(result):
            handle.write("\t".join(fmt(row.get(c)) for c in FLAG_COLUMNS) + "\n")
    keys = (
        "mark", "tier", "samples", "groups", "singletons", "flagged", "median_r_within_group",
        "median_r_between_groups", "batch_delta", "low_quality_floor", "note",
    )  # fmt: skip
    with open(destination / "qc_marks.tsv", "w", encoding="utf-8") as handle:
        handle.write("\t".join(keys) + "\n")
        for info in result.mark_summary:
            handle.write("\t".join(fmt(info.get(k)) for k in keys) + "\n")
    labels = sorted(result.confusion)
    predicted = sorted({p for c in result.confusion.values() for p in c})
    with open(destination / "qc_mark_confusion.tsv", "w", encoding="utf-8") as handle:
        handle.write("labelled\\predicted\t" + "\t".join(predicted) + "\n")
        for m in labels:
            handle.write(
                m + "\t" + "\t".join(str(result.confusion[m].get(p, 0)) for p in predicted) + "\n"
            )
    with open(destination / "qc_group_matrix.tsv", "w", encoding="utf-8") as handle:
        handle.write("mark\tgroup_a\tgroup_b\tmedian_r\tpairs\n")
        for mark, (names, medians, pairs) in result.group_matrices.items():
            for p, g in enumerate(names):
                for q, h in enumerate(names):
                    value = "" if np.isnan(medians[p, q]) else f"{medians[p, q]:.4f}"
                    handle.write(f"{mark}\t{g}\t{h}\t{value}\t{pairs[p, q]}\n")


def qc_from_files(sheet: str | Path, matrix: str | Path, options: dict) -> tuple[QCResult, dict]:
    """Map sheet columns, apply filters, match rows to matrix columns and run :func:`run_qc`.

    Options (all optional except the columns the guesses cannot fill): ``qc_sample``,
    ``qc_mark``, ``qc_group``, ``qc_batch``, ``qc_filters``, ``qc_k``, ``qc_top``,
    ``qc_transform``, ``qc_method``, ``qc_z``, ``qc_margin``, ``qc_skip_marks``.
    """
    columns, rows = read_sheet(sheet)
    data = load_matrix(matrix)
    guesses = guess_columns(columns, rows, data.samples)
    chosen = {
        k: options.get(f"qc_{k}") or guesses.get(k) for k in ("sample", "mark", "group", "batch")
    }
    for role in ("sample", "mark", "group"):
        if not chosen[role]:
            raise ValueError(f"Choose the sample-sheet column that holds the {role}")
    for role, column in chosen.items():
        if column and column not in columns:
            raise ValueError(f"Column {column!r} ({role}) is not in the sample sheet")
    kept = apply_filters(rows, parse_filters(options.get("qc_filters", "")))
    matched = match_samples(kept, chosen["sample"], data.samples)
    if len(matched) < 3:
        raise ValueError(
            f"Only {len(matched)} of {len(kept)} sheet rows match matrix columns by "
            f"{chosen['sample']!r}; pick the column whose values equal the matrix sample names"
        )
    order = sorted(matched, key=matched.get)
    columns_used = [matched[i] for i in order]
    if len(set(columns_used)) != len(columns_used):
        raise ValueError(
            f"Several sheet rows map to the same matrix column via {chosen['sample']!r}"
        )
    meta = [
        dict(
            sample=kept[i][chosen["sample"]],
            mark=_value(kept[i].get(chosen["mark"])),
            group=_value(kept[i].get(chosen["group"])),
            batch=_value(kept[i].get(chosen["batch"])) if chosen["batch"] else "",
        )
        for i in order
    ]
    subset = SignalMatrix(
        [data.samples[j] for j in columns_used],
        data.values[:, columns_used],
        data.chroms,
        data.starts,
        data.bin_size,
    )
    skip = tuple(
        s.strip() for s in str(options.get("qc_skip_marks", "INPUT, IgG")).split(",") if s.strip()
    )
    result = run_qc(
        subset,
        meta,
        k=int(options.get("qc_k") or 5),
        top=int(options.get("qc_top") or 5000),
        transform=str(options.get("qc_transform") or "log1p"),
        method=str(options.get("qc_method") or "pearson"),
        z=float(options.get("qc_z") or 3.0),
        margin=float(
            options.get("qc_margin") if options.get("qc_margin") not in (None, "") else 0.02
        ),
        skip_group_marks=skip,
    )
    unmatched = [kept[i].get(chosen["sample"], "") for i in range(len(kept)) if i not in matched]
    notes = dict(
        columns=chosen,
        sheet_rows=len(rows),
        after_filters=len(kept),
        matched=len(matched),
        unmatched=unmatched[:50],
        matrix_samples=len(data.samples),
        features=int(data.values.shape[0]),
        bin_size=data.bin_size,
        skip_group_marks=list(skip),
    )
    return result, notes


def _value(raw) -> str:
    text = str(raw or "").strip()
    return "" if text.lower() in MISSING else text
