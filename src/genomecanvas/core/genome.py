"""Genome-level coordinate metadata."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Genome:
    """A named genome assembly with optional chromosome lengths."""

    assembly: str
    chromosomes: tuple[tuple[str, int], ...] = ()
