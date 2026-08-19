"""Materialized visualization-derived signal tracks."""

from __future__ import annotations

from dataclasses import dataclass

from ..model.background_spec import BackgroundSpec
from ..model.track import SignalSlice, Track
from ..transform.background import (
    BackgroundEstimate,
    estimate_background,
    transform_values,
)


@dataclass(frozen=True, slots=True)
class DerivedTrack:
    """A derived signal with its parent, transformation, and audit metadata."""

    id: str
    label: str
    parent: Track
    background_spec: BackgroundSpec
    operation: str
    pseudocount: float | None
    background: BackgroundEstimate
    result: SignalSlice
    parents: tuple[str, ...] = ()

    @property
    def type(self) -> str:
        """Return the declarative track type."""
        return "derived_signal"

    @property
    def metadata(self) -> dict[str, object]:
        """Return provenance-ready transformation metadata."""
        return {
            "parent_track_id": self.parent.id,
            "parent_track_ids": self.parents or (self.parent.id,),
            "background_strategy": self.background_spec.strategy,
            "background_estimator": self.background.estimator,
            "upstream": self.background_spec.upstream,
            "downstream": self.background_spec.downstream,
            "background_value": self.background.value,
            "left_background_value": self.background.left_value,
            "right_background_value": self.background.right_value,
            "transformation": self.operation,
            "pseudocount": self.pseudocount,
            "excluded_intervals": self.background_spec.exclusions,
        }


def derive_signal(
    parent: Track,
    target: SignalSlice,
    left_background: SignalSlice,
    right_background: SignalSlice,
    background_spec: BackgroundSpec,
    operation: str,
    pseudocount: float | None = None,
) -> DerivedTrack:
    """Materialize a background-adjusted signal track once."""
    estimate = estimate_background(
        left_background.values, right_background.values, background_spec
    )
    values = transform_values(target.values, estimate.value, operation, pseudocount)
    result = SignalSlice(
        target.chrom, target.start, target.end, target.bin_starts, target.bin_ends,
        values, target.units, {"derived": "local_background"},
    )
    return DerivedTrack(
        f"{parent.id}_{operation}_local_background", f"{parent.label} · local {operation}",
        parent, background_spec, operation, pseudocount, estimate, result, (parent.id,),
    )
