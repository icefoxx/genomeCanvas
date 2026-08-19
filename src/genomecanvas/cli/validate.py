"""Validation command implementation."""

from ..io.validation import validate_scientific
from ..model.figure_spec import FigureSpec


def run(path: str) -> list[str]:
    issues = validate_scientific(FigureSpec.from_yaml(path))
    return [f"{issue.severity}: {issue.message}" for issue in issues]
