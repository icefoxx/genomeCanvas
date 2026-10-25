"""Declarative overlays and click-edit payloads reserved for interactive UIs."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field, validator

OverlayKind = Literal[
    "label", "text", "highlight", "arrow", "region_note", "point_note"
]
EditAction = Literal[
    "add_overlay",
    "update_overlay",
    "remove_overlay",
    "update_track_style",
    "set_track_height",
    "set_track_visibility",
    "update_figure_style",
    "set_title",
]


class OverlayAnnotation(BaseModel):
    """A user-facing note drawn on the figure, addressable by id for later edits."""

    id: str
    kind: OverlayKind = "label"
    chrom: str | None = None
    start: int | None = None
    end: int | None = None
    x: float | None = None
    y: float | None = None
    text: str | None = None
    track_id: str | None = None
    color: str = "#111111"
    fill: str = "#ffc41a"
    alpha: float = Field(default=0.35, ge=0, le=1)
    font_size: float | None = None
    font_style: Literal["normal", "italic", "oblique"] | None = None
    font_weight: Literal["normal", "bold"] | None = None
    ha: Literal["left", "center", "right"] = "left"
    va: Literal["top", "center", "bottom"] = "bottom"

    @validator("id")
    def nonempty(cls, value: str) -> str:
        if not value:
            raise ValueError("overlay id must not be empty")
        return value


class FigureEdit(BaseModel):
    """One atomic mutation a click UI can apply to a FigureSpec."""

    action: EditAction
    overlay: OverlayAnnotation | None = None
    overlay_id: str | None = None
    track_id: str | None = None
    style: dict[str, Any] | None = None
    height: float | None = None
    visibility: bool | None = None
    title: str | None = None


class HitTarget(BaseModel):
    """A pickable object recorded at render time for future click mapping."""

    id: str
    kind: str
    track_id: str | None = None
    overlay_id: str | None = None
    chrom: str
    start: int
    end: int
    text: str | None = None
    axes_index: int
    x: float
    y: float
    payload: dict[str, Any] = Field(default_factory=dict)
