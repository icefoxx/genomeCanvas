"""Typed figure and font styling for locus plots."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

FontStyle = Literal["normal", "italic", "oblique"]
FontWeight = Literal["normal", "bold"]


class FontSpec(BaseModel):
    """A reusable font declaration."""

    family: str | None = None
    size: float = Field(default=8.0, gt=0)
    color: str = "#222222"
    style: FontStyle = "normal"
    weight: FontWeight = "normal"

    def matplotlib(self, family: str | None = None) -> dict[str, Any]:
        """Return kwargs for ``Axes.text`` / ``set_ylabel`` / ``tick_params``."""
        return {
            "fontsize": self.size,
            "color": self.color,
            "fontstyle": self.style,
            "fontweight": self.weight,
            "fontfamily": self.family or family or "DejaVu Sans",
        }


class FigureStyle(BaseModel):
    """Figure-wide typography and chrome."""

    font_family: str = "DejaVu Sans"
    title: str | None = None
    title_font: FontSpec = Field(
        default_factory=lambda: FontSpec(size=11, weight="bold")
    )
    axis_label_font: FontSpec = Field(default_factory=lambda: FontSpec(size=9))
    tick_font: FontSpec = Field(default_factory=lambda: FontSpec(size=7))
    track_label_font: FontSpec = Field(default_factory=lambda: FontSpec(size=8))
    gene_label_font: FontSpec = Field(
        default_factory=lambda: FontSpec(size=7, style="italic", color="#444444")
    )
    legend_font: FontSpec = Field(default_factory=lambda: FontSpec(size=7))
    overlay_font: FontSpec = Field(
        default_factory=lambda: FontSpec(size=8, weight="bold", color="#111111")
    )
    background: str = "white"
    spine_color: str = "#444444"
    grid: bool = False
    signal_fill_alpha: float = Field(default=0.85, ge=0, le=1)
    signal_linewidth: float = Field(default=0.6, ge=0)


DEFAULT_TRACK_HEIGHTS: dict[str, float] = {
    "signal": 1.0,
    "derived_signal": 1.0,
    "difference": 1.0,
    "overlay": 1.3,
    "peak": 0.4,
    "peak_signal": 0.4,
    "gene": 0.9,
    "isoform": 1.2,
    "chromhmm": 0.35,
    "cytoband": 0.4,
    "loop": 0.85,
    "variant": 0.65,
    "motif": 0.5,
    "annotation": 0.45,
    "heatmap": 1.4,
    "highlight": 0.25,
}


def track_height(track_type: str, explicit: float | None) -> float:
    """Resolve a track height from an explicit value or a type default."""
    if explicit is not None and explicit > 0:
        return explicit
    return DEFAULT_TRACK_HEIGHTS.get(track_type, 1.0)


def font_from_style(style: dict[str, Any], base: FontSpec, family: str) -> FontSpec:
    """Overlay optional per-track font keys onto a figure-level font."""
    return FontSpec(
        family=style.get("font_family", base.family or family),
        size=float(style.get("font_size", base.size)),
        color=str(style.get("font_color", base.color)),
        style=style.get("font_style", base.style),
        weight=style.get("font_weight", base.weight),
    )
