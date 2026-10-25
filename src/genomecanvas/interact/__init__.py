"""Reserved interactive-edit API: overlays, hit targets, and FigureEdit application."""

from ..model.overlays import FigureEdit, HitTarget, OverlayAnnotation
from .edits import apply_edit, overlay_by_id
from .hits import HitRegistry, hit_at_genome

__all__ = [
    "FigureEdit",
    "HitRegistry",
    "HitTarget",
    "OverlayAnnotation",
    "apply_edit",
    "hit_at_genome",
    "overlay_by_id",
]
