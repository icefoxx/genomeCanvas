"""Scientific figure-choice linting for advanced transformations."""

from __future__ import annotations

from ..io.validation import ValidationIssue
from ..model.figure_spec import FigureSpec


def lint_figure(spec: FigureSpec) -> list[ValidationIssue]:
    """Return errors, warnings, and informational diagnostics."""
    issues: list[ValidationIssue] = []
    signal_tracks = [track for track in spec.tracks if track.type == "signal"]
    if len(signal_tracks) > 1 and spec.scales.get("mode") == "per_track":
        issues.append(ValidationIssue("WARNING", "independent y-axes may hinder signal comparison"))
    for transformation in spec.transformations:
        if (
            transformation.operation in {"ratio", "log2_ratio"}
            and transformation.pseudocount is not None
        ):
            if transformation.pseudocount > 1:
                issues.append(
                    ValidationIssue(
                        "WARNING", f"{transformation.id}: pseudocount may dominate signal"
                    )
                )
        if transformation.background is not None and any(
            track.metadata.get("broad_mark") for track in signal_tracks
        ):
            issues.append(
                ValidationIssue(
                    "WARNING", f"{transformation.id}: broad marks may invalidate local background"
                )
            )
    for track in signal_tracks:
        if track.metadata.get("empty"):
            issues.append(ValidationIssue("WARNING", f"track {track.id!r} is empty"))
    return issues
