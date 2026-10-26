"""Inputs for ProfileSeg: config, sample manifest, regions, and the profile tensor."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

MANIFEST_COLUMNS = [
    "sample_id",
    "condition",
    "biological_replicate",
    "assay",
    "genome_build",
    "signal_file",
    "signal_type",
    "normalization",
    "region_set",
    "batch",
    "pairing_id",
    "bigwig_params",
    "fragment_length",
]
ANCHORS = ("center", "start", "end", "summit")


# ---------------------------------------------------------------- configuration


@dataclass
class Config:
    manifest: str = "metadata/sample_manifest.tsv"
    regions: str = ""
    matrix: str = (
        ""  # precomputed .npz (ProfileSeg) or deepTools computeMatrix .gz instead of bigWigs
    )
    reference: str = "WT"
    alternative: str = "KO"
    anchor: str = "center"
    orientation: str = "auto"  # auto | oriented | unoriented
    flank: int = 2000
    binsize: int = 20
    background_width: int = 500
    background_estimator: str = (
        "mean"  # mean (unbiased for count noise) | median (robust to flank peaks)
    )
    core_flank: int = 1000
    negative_policy: str = "signed"  # signed (moments on the signed aggregate) | exclude
    max_truncated_fraction: float = 0.2
    min_mass_quantile: float = 0.1
    # shrunk_sample (default): per-sample bounded scaling (mass + sample median); shrunk: pooled
    # across samples; equal: pooled mass only; signal: no per-region scaling.
    region_weighting: str = "shrunk_sample"
    # Condition-blind shape classes (shift-only k-mean alignment of pooled profiles, after FunChIP):
    # an optional descriptive table of which kinds of profiles change. Classes do not enter any test
    # (stratifying the pooled nulls by class was evaluated and rejected: no calibration gain, less power).
    shape_classes: bool = False
    # Re-centre each region on its class template before all analyses (corrects anchor misplacement;
    # the shift comes from the pooled profile and is applied to every sample alike). Off by default;
    # recommended for imprecise anchors (summits, motifs) and aggregate statements.
    align_anchors: bool = False
    # Manifest `batch` column (when filled for every sample). batch_aware: refuse per-region shape tests
    # when condition and batch coincide (no assignment within batches exists). batch_permutation:
    # "free" (default) keeps all label assignments for balanced designs; "within" permutes labels only
    # within batches (more power for balanced designs, FDP ~0.07 vs ~0.03 on real yeast data).
    batch_aware: bool = True
    batch_permutation: str = "free"
    align_max_shift: int = 150
    center_half: int = 250
    primary_amplitude: str = "log2_auc"
    primary_shape: str = "auto"  # auto: centroid for oriented anchors, iqr_width otherwise
    bootstrap: int = 1000
    region_fdr: float = 0.05
    confidence: float = 0.95
    seed: int = 7201
    output: str = "ProfileSeg_out"
    palette: list = field(default_factory=lambda: ["#0072B2", "#D55E00"])
    font: str = "DejaVu Sans"
    dpi: int = 300
    base_dir: str = "."

    def path(self, value: str) -> Path:
        p = Path(value).expanduser()
        return p if p.is_absolute() else Path(self.base_dir) / p


def load_config(path: str | Path, overrides: dict | None = None) -> Config:
    path = Path(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    flat: dict = {}
    for key, value in data.items():
        if isinstance(value, dict):  # allow grouped sections: contrast:, background:, shape:, ...
            for sub, item in value.items():
                flat[_ALIASES.get(f"{key}.{sub}", sub)] = item
        else:
            flat[key] = value
    flat.update({k: v for k, v in (overrides or {}).items() if v is not None})
    known = Config.__dataclass_fields__
    unknown = sorted(set(flat) - set(known))
    if unknown:
        raise ValueError(f"Unknown ProfileSeg setting(s): {', '.join(unknown)}")
    config = Config(**flat)
    config.base_dir = str(path.parent.resolve())
    validate_config(config)
    return config


_ALIASES = {
    "contrast.reference": "reference",
    "contrast.alternative": "alternative",
    "background.width": "background_width",
    "background.estimator": "background_estimator",
    "shape.core_flank": "core_flank",
    "shape.negative_policy": "negative_policy",
    "shape.max_truncated_fraction": "max_truncated_fraction",
    "shape.min_mass_quantile": "min_mass_quantile",
    "shape.region_weighting": "region_weighting",
    "shape.classes": "shape_classes",
    "shape.align_anchors": "align_anchors",
    "shape.align_max_shift": "align_max_shift",
    "inference.batch_aware": "batch_aware",
    "inference.batch_permutation": "batch_permutation",
    "amplitude.center_half": "center_half",
    "inference.primary_amplitude": "primary_amplitude", "inference.primary_shape": "primary_shape",
    "inference.bootstrap": "bootstrap", "inference.region_fdr": "region_fdr",
    "inference.confidence": "confidence",
    "inference.seed": "seed",
    "figures.palette": "palette",
    "figures.font": "font",
    "figures.dpi": "dpi",
}


def validate_config(c: Config) -> None:
    if c.anchor not in ANCHORS:
        raise ValueError(f"anchor must be one of {ANCHORS}")
    if c.orientation not in ("auto", "oriented", "unoriented"):
        raise ValueError("orientation must be auto, oriented or unoriented")
    if c.primary_amplitude not in ("log2_auc", "log2_centre_auc", "log2_enrichment"):
        raise ValueError("primary_amplitude must be log2_auc, log2_centre_auc or log2_enrichment")
    if c.primary_shape not in ("auto", "centroid", "sd_width", "iqr_width", "asymmetry"):
        raise ValueError("primary_shape must be auto, centroid, sd_width, iqr_width or asymmetry")
    if c.region_weighting not in ("shrunk", "shrunk_sample", "equal", "signal"):
        raise ValueError("region_weighting must be shrunk, shrunk_sample, equal or signal")
    if not 0 < c.align_max_shift < c.core_flank:
        raise ValueError("shape.align_max_shift must be positive and smaller than core_flank")
    if c.batch_permutation not in ("free", "within"):
        raise ValueError("inference.batch_permutation must be free or within")
    if c.bootstrap < 10:
        raise ValueError("inference.bootstrap must be at least 10 (bootstrap intervals are reported for every metric)")
    if c.background_estimator not in ("mean", "median"):
        raise ValueError("background estimator must be mean or median")
    if c.negative_policy == "truncate":  # v1 drafts used this name
        c.negative_policy = "signed"
    if c.negative_policy not in ("signed", "exclude"):
        raise ValueError("negative_policy must be signed or exclude")
    if c.flank % c.binsize:
        raise ValueError("flank must be a multiple of binsize")
    if not 0 < c.background_width <= c.flank or not 0 < c.core_flank <= c.flank:
        raise ValueError("background_width and core_flank must be within (0, flank]")
    if c.core_flank > c.flank - c.background_width:
        raise ValueError(
            "core window overlaps the background flanks; reduce core_flank or background_width"
        )
    if c.reference == c.alternative:
        raise ValueError("reference and alternative conditions must differ")


EXAMPLE_CONFIG = """# ProfileSeg configuration (v1)
manifest: metadata/sample_manifest.tsv   # columns: see docs/PROFILESEG_PLAN.md
regions: regions.bed                     # fixed, condition-independent regions (BED6; strand for oriented anchors)
# matrix: precomputed.npz                # alternative input instead of bigWigs + regions
contrast: {reference: WT, alternative: KO}
anchor: center                           # center | start (TSS) | end (TES) | summit (narrowPeak col 10)
orientation: auto                        # auto: oriented for start/end anchors with +/- strands
flank: 2000
binsize: 20
background: {width: 500, estimator: mean}  # outer flanks; mean is unbiased for count noise
shape:
  core_flank: 1000                       # shape metrics use ±core_flank after background subtraction
  negative_policy: signed                # signed aggregate moments (no clipping) | exclude
  max_truncated_fraction: 0.2
  min_mass_quantile: 0.1                 # condition-independent: pooled mass across all samples
  region_weighting: shrunk_sample        # per-sample scale = mass + sample median (bounded, label-free)
amplitude: {center_half: 250}
inference:
  primary_amplitude: log2_auc            # pre-specified; other metrics are Holm-adjusted
  primary_shape: auto                    # auto: centroid (oriented) | iqr_width (unoriented)
  bootstrap: 1000
  confidence: 0.95
  seed: 7201
output: ProfileSeg_out
figures: {palette: ["#0072B2", "#D55E00"], font: DejaVu Sans, dpi: 300}
"""


# ---------------------------------------------------------------- manifest


def read_manifest(path: Path) -> tuple[list[dict], list[str]]:
    """Rows plus human-readable flags (missing metadata is flagged, never guessed)."""
    rows = []
    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        rows.append({k.strip(): (v or "").strip() for k, v in row.items() if k})
    flags: list[str] = []
    if not rows:
        raise ValueError(f"{path} has no samples")
    missing_cols = [c for c in MANIFEST_COLUMNS if c not in rows[0]]
    if missing_cols:
        flags.append("manifest lacks columns: " + ", ".join(missing_cols))
    for required in ("sample_id", "condition"):
        if any(not r.get(required) for r in rows):
            raise ValueError(f"every manifest row needs '{required}'")
    if len({r["sample_id"] for r in rows}) != len(rows):
        raise ValueError("sample_id values must be unique")
    builds = {r.get("genome_build", "") for r in rows} - {""}
    if len(builds) > 1:
        raise ValueError(f"mixed genome builds in manifest: {sorted(builds)}")
    if not builds:
        flags.append("genome_build is not recorded")
    for column in (
        "biological_replicate",
        "normalization",
        "bigwig_params",
        "fragment_length",
        "assay",
    ):
        empty = [r["sample_id"] for r in rows if not r.get(column)]
        if empty:
            flags.append(f"'{column}' missing for {len(empty)} sample(s): {', '.join(empty[:6])}")
    fragments: dict[str, list[float]] = {}
    for r in rows:
        try:
            value = float(r.get("fragment_length", ""))
        except ValueError:
            continue
        fragments.setdefault(r["condition"], []).append(value)
    if len(fragments) >= 2:
        med = {k: float(np.median(v)) for k, v in fragments.items()}
        if max(med.values()) > 1.2 * min(med.values()):
            flags.append(
                f"fragment length differs between conditions {med}: apparent width changes may be technical"
            )
    return rows, flags


# ---------------------------------------------------------------- profiles


@dataclass
class Profiles:
    """Signal tensor (samples × regions × bins) on physical bp edges relative to the anchor."""

    values: np.ndarray
    edges: np.ndarray
    samples: list[dict]
    regions: list[dict]
    oriented: bool
    provenance: dict = field(default_factory=dict)

    @property
    def conditions(self) -> np.ndarray:
        return np.array([s["condition"] for s in self.samples])

    def save(self, path: Path) -> None:
        np.savez_compressed(
            path,
            values=self.values,
            edges=self.edges,
            oriented=self.oriented,
            samples=json.dumps(self.samples),
            regions=json.dumps(self.regions),
            provenance=json.dumps(self.provenance),
        )

    @classmethod
    def load(cls, path: Path) -> Profiles:
        with np.load(path, allow_pickle=False) as data:
            return cls(
                values=data["values"].astype(float),
                edges=data["edges"].astype(float),
                samples=json.loads(str(data["samples"])),
                regions=json.loads(str(data["regions"])),
                oriented=bool(data["oriented"]),
                provenance=json.loads(str(data["provenance"])),
            )


def read_regions(path: Path) -> list[dict]:
    regions = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                raise ValueError(f"BED line needs ≥3 columns: {line.strip()}")
            regions.append(
                dict(
                    chrom=f[0],
                    start=int(f[1]),
                    end=int(f[2]),
                    name=f[3] if len(f) > 3 and f[3] else f"r{len(regions) + 1}",
                    strand=f[5] if len(f) > 5 and f[5] in "+-" and f[5] else ".",
                    summit=int(f[9])
                    if len(f) > 9 and f[9].lstrip("-").isdigit() and int(f[9]) >= 0
                    else None,
                )
            )
    if not regions:
        raise ValueError(f"{path} contains no regions")
    return regions


def anchor_base(region: dict, anchor: str) -> int:
    """0-based anchor base; 'center' rounds toward the 5' end so ± windows are exact mirrors."""
    s, e, strand = region["start"], region["end"], region["strand"]
    if anchor == "summit":
        if region.get("summit") is None:
            raise ValueError(f"region {region['name']} has no summit offset (narrowPeak column 10)")
        return s + region["summit"]
    if anchor == "start":
        return e - 1 if strand == "-" else s
    if anchor == "end":
        return s if strand == "-" else e - 1
    return (s + e - 1) // 2 if strand == "-" else (s + e) // 2


def resolve_orientation(config: Config, regions: list[dict]) -> bool:
    if config.orientation == "oriented":
        if any(r["strand"] not in "+-" or r["strand"] == "." for r in regions):
            raise ValueError("orientation: oriented requires a +/- strand for every region")
        return True
    if config.orientation == "unoriented":
        return False
    stranded = all(r["strand"] in ("+", "-") for r in regions)
    return stranded and config.anchor in ("start", "end")


def extract(config: Config, samples: list[dict], regions: list[dict]) -> Profiles:
    """bigWig/bedGraph → tensor using the GenomeCanvas window reader (NaN beyond chromosome ends)."""
    from ..pipeline.genome_bins import window_matrix

    oriented = resolve_orientation(config, regions)
    flank, binsize = config.flank, config.binsize
    n_bins = 2 * flank // binsize
    windows, minus = [], []
    for region in regions:
        base = anchor_base(region, config.anchor)
        neg = oriented and region["strand"] == "-"
        lo = base + 1 - flank if neg else base - flank  # mirror around the anchor base on '-'
        windows.append((region["chrom"], lo, lo + 2 * flank))
        minus.append(neg)
    minus = np.array(minus)
    values = np.empty((len(samples), len(regions), n_bins))
    for i, sample in enumerate(samples):
        path = config.path(sample["signal_file"])
        if not path.is_file():
            raise FileNotFoundError(f"signal file for {sample['sample_id']} not found: {path}")
        block = window_matrix(str(path), windows, n_bins, "mean")
        block[minus] = block[minus, ::-1]
        values[i] = block
    edges = np.arange(-flank, flank + binsize, binsize, dtype=float)
    return Profiles(
        values,
        edges,
        samples,
        regions,
        oriented,
        provenance=dict(source="bigwig", anchor=config.anchor, flank=flank, binsize=binsize),
    )


def read_deeptools_matrix(path: Path, samples: list[dict]) -> Profiles:
    """deepTools computeMatrix (reference-point) output: header line '@{json}' then BED6 + values."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        header = json.loads(handle.readline().lstrip("@"))
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    labels = header["sample_labels"]
    bounds = header["sample_boundaries"]
    up = int(np.atleast_1d(header["upstream"])[0])
    down = int(np.atleast_1d(header["downstream"])[0])
    binsize = int(np.atleast_1d(header["bin size"])[0])
    if int(np.atleast_1d(header.get("body", [0]))[0]):
        raise ValueError("scale-regions matrices are not supported in v1; use reference-point")
    data = np.array(
        [[float(x) if x not in ("nan", "NA", "") else np.nan for x in r[6:]] for r in rows]
    )
    by_label = {s.get("matrix_label") or s["sample_id"]: s for s in samples}
    values, chosen = [], []
    for k, label in enumerate(labels):
        if label in by_label:
            values.append(data[:, bounds[k] : bounds[k + 1]])
            chosen.append(by_label[label])
    if len(chosen) != len(samples):
        raise ValueError(f"matrix labels {labels} do not match manifest sample_id/matrix_label")
    regions = [
        dict(
            chrom=r[0],
            start=int(r[1]),
            end=int(r[2]),
            name=r[3],
            strand=r[5] if r[5] in "+-" else ".",
        )
        for r in rows
    ]
    edges = np.arange(-up, down + binsize, binsize, dtype=float)
    # deepTools already orients minus-strand rows.
    oriented = all(r["strand"] in "+-" for r in regions) and header.get("ref point", ["center"])[
        0
    ] in ("TSS", "TES")
    return Profiles(
        np.stack(values),
        edges,
        chosen,
        regions,
        oriented,
        provenance=dict(source="deeptools", file=str(path)),
    )
