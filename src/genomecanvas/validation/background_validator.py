"""Validation and lint rules for local-background transformations."""

from __future__ import annotations

from dataclasses import dataclass

from ..model.background_spec import BackgroundSpec
from ..transform.background import BackgroundEstimate


@dataclass(frozen=True, slots=True)
class BackgroundIssue:
    """Severity-tagged background diagnostic."""

    severity: str
    message: str


def validate_background(
    spec: BackgroundSpec, estimate: BackgroundEstimate | None = None
) -> list[BackgroundIssue]:
    """Validate explicit flanks and report local-background risks."""
    issues: list[BackgroundIssue] = []
    if (
        (spec.upstream is not None and spec.upstream < 0)
        or (spec.downstream is not None and spec.downstream < 0)
    ):
        issues.append(BackgroundIssue("ERROR", "flank lengths cannot be negative"))
    if spec.strategy == "flanking" and spec.upstream == 0 and spec.downstream == 0:
        issues.append(
            BackgroundIssue("ERROR", "background is unavailable because both flanks are empty")
        )
    if len(spec.exclusions) > 10:
        issues.append(BackgroundIssue("WARNING", "background contains many excluded intervals"))
    if estimate is not None:
        issues.extend(BackgroundIssue("WARNING", message) for message in estimate.warnings)
        if abs(estimate.value) < 1e-12:
            issues.append(
                BackgroundIssue("WARNING", "background close to zero makes ratios unstable")
            )
    return issues
