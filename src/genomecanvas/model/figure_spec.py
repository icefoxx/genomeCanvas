"""Typed, serializable declarative figure specification."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, Field, root_validator, validator

from .background_spec import TransformationSpec


class GenomeSpec(BaseModel):
    """Genome assembly declaration."""

    assembly: str


class RegionSpec(BaseModel):
    """0-based, half-open figure region."""

    chrom: str
    start: int
    end: int

    @validator("chrom")
    def nonempty_chrom(cls, value: str) -> str:
        if not value:
            raise ValueError("region chrom must not be empty")
        return value

    @root_validator(skip_on_failure=True)
    def valid_bounds(cls, values: dict[str, Any]) -> dict[str, Any]:
        if values.get("start", -1) < 0 or values.get("end", 0) <= values.get("start", 0):
            raise ValueError("region must satisfy 0 <= start < end")
        return values


class SourceSpec(BaseModel):
    """Input source declaration."""

    path: str


class TrackSpec(BaseModel):
    """Declarative track metadata."""

    id: str
    type: str
    label: str | None = None
    source: SourceSpec | None = None
    height: float = Field(default=1.0, gt=0)
    visibility: bool = True
    metadata: dict[str, Any] = Field(default_factory=dict)
    transformations: list[TransformationSpec] = Field(default_factory=list)
    assay: str | None = None
    condition: str | None = None
    replicate: str | None = None
    sample: str | None = None
    group: str | None = None
    units: str | None = None
    genome_assembly: str | None = None
    scale_group: str | None = None
    z_order: int = 0
    style: dict[str, Any] = Field(default_factory=dict)

    @validator("id")
    def nonempty_id(cls, value: str) -> str:
        if not value:
            raise ValueError("track id must not be empty")
        return value


class LayoutSpec(BaseModel):
    """Figure dimensions."""

    width: float = Field(default=7.0, gt=0)
    height: float | None = Field(default=None, gt=0)
    aspect_ratio: float | None = Field(default=None, gt=0)
    track_gap: float = Field(default=0.02, ge=0, lt=1)
    upper_annotation_ratio: float | None = Field(default=None, gt=0)
    units: Literal["inch", "mm"] = "inch"
    orientation: Literal["vertical", "horizontal"] = "vertical"
    x_label: str | None = None
    show_x_axis: bool = True
    show_gene_track: bool = True
    show_cytoband: bool = False
    cytoband_marker: Literal["bar", "line", "triangle", "circle", "diamond"] = "bar"
    gene_track_scale: float = Field(default=1.0, gt=0)
    cytoband_scale: float = Field(default=1.0, gt=0)
    gene_cytoband_gap: float = Field(default=0.02, ge=0)
    cytoband_marker_color: str = "#e33d32"
    cytoband_marker_shape: Literal["bar", "line", "triangle", "circle", "diamond"] = "bar"
    gene_model_label: Literal["none", "name", "id"] = "none"
    transcript_label: Literal["id", "name", "none"] = "id"
    x_tick_count: int = Field(default=5, ge=2, le=20)
    y_axis_mode: Literal["outer", "all", "left", "right"] = "outer"
    y_value_transform: Literal["identity", "divide", "log2", "log10"] = "identity"
    y_divisor: float = Field(default=1.0, gt=0)
    y_tick_count: int = Field(default=5, ge=2, le=20)
    y_tick_decimals: int = Field(default=1, ge=0, le=8)
    y_tick_label_pad: float = 3.0
    show_y_ticks: bool = True
    show_y_axis: bool = True
    y_tick_length: float = Field(default=3.0, ge=0)
    show_major_frame: bool = True


def _default_formats() -> list[Literal["pdf", "svg", "png"]]:
    """Return default vector export formats."""
    return ["pdf", "svg"]


class ExportSpec(BaseModel):
    """Requested output formats."""

    formats: list[Literal["pdf", "svg", "png"]] = Field(default_factory=_default_formats)


class FigureSpec(BaseModel):
    """Top-level validated figure declaration."""

    version: str = "1.0"
    genome: GenomeSpec
    region: RegionSpec
    tracks: list[TrackSpec] = Field(default_factory=list)
    groups: list[dict[str, Any]] = Field(default_factory=list)
    scales: dict[str, Any] = Field(default_factory=dict)
    highlights: list[dict[str, Any]] = Field(default_factory=list)
    zoom: dict[str, Any] | None = None
    layout: LayoutSpec = Field(default_factory=LayoutSpec)
    theme: str = "manuscript"
    export: ExportSpec = Field(default_factory=ExportSpec)
    metadata: dict[str, Any] = Field(default_factory=dict)
    transformations: list[TransformationSpec] = Field(default_factory=list)

    @root_validator(skip_on_failure=True)
    def unique_track_ids(cls, values: dict[str, Any]) -> dict[str, Any]:
        tracks = values.get("tracks", [])
        ids = [track.id for track in tracks]
        if len(ids) != len(set(ids)):
            raise ValueError("track IDs must be unique")
        return values

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON/YAML-compatible representation."""
        return self.dict(exclude_none=True)

    def to_yaml(self, path: str | Path | None = None) -> str:
        """Serialize this specification to YAML, optionally writing a file."""
        text = yaml.safe_dump(self.to_dict(), sort_keys=False)
        if path is not None:
            Path(path).write_text(text, encoding="utf-8")
        return text

    @classmethod
    def from_yaml(cls, path: str | Path) -> FigureSpec:
        """Load and validate YAML or JSON from a file."""
        source_path = Path(path).resolve()
        payload = yaml.safe_load(source_path.read_text(encoding="utf-8"))
        for track in payload.get("tracks", []):
            source = track.get("source")
            if source and "path" in source:
                source["path"] = str((source_path.parent / source["path"]).resolve())
        return cls.parse_obj(payload)

    def to_json(self) -> str:
        """Serialize this specification to stable JSON."""
        return json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))
