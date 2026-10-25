"""Schema-adjacent scientific validation diagnostics."""

from dataclasses import dataclass

from ..model.figure_spec import FigureSpec
from ..validation.background_validator import validate_background


@dataclass(frozen=True, slots=True)
class ValidationIssue:
    severity: str
    message: str


def validate_scientific(spec: FigureSpec) -> list[ValidationIssue]:
    """Check explicit assembly and comparison assumptions."""
    issues: list[ValidationIssue] = []
    assemblies = {track.genome_assembly for track in spec.tracks if track.genome_assembly}
    if len(assemblies) > 1:
        issues.append(ValidationIssue("ERROR", "tracks declare incompatible genome assemblies"))
    for track in spec.tracks:
        if track.type == "overlay" and not track.members:
            issues.append(
                ValidationIssue("ERROR", f"track {track.id!r} overlay is missing members")
            )
        sourced = {
            "signal", "peak", "peak_signal", "gene", "isoform", "chromhmm", "cytoband",
            "annotation", "loop", "variant", "motif",
        }
        gene_ok = track.type in {"gene", "isoform"} and spec.gene_model
        if track.type in sourced and track.source is None and not gene_ok:
            issues.append(ValidationIssue("ERROR", f"track {track.id!r} is missing a source"))
    for transformation in spec.transformations:
        if transformation.background is not None:
            for issue in validate_background(transformation.background):
                issues.append(
                    ValidationIssue(issue.severity, f"{transformation.id}: {issue.message}")
                )
        if (
            transformation.operation in {"ratio", "log2_ratio"}
            and transformation.pseudocount is None
        ):
            issues.append(
                ValidationIssue(
                    "ERROR",
                    f"{transformation.id}: ratio operations require an explicit pseudocount",
                )
            )
    return issues
