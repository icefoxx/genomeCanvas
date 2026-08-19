"""Explicit source-to-internal coordinate conversions."""

from __future__ import annotations

from dataclasses import dataclass

from .errors import CoordinateError

INTERNAL_COORDINATE_SYSTEM = "0-based, half-open"


@dataclass(frozen=True, slots=True)
class CoordinateConversion:
    """Document an explicit conversion applied by an input adapter."""

    source_coordinate_system: str
    internal_coordinate_system: str = INTERNAL_COORDINATE_SYSTEM
    conversion_applied: str = "none"


def one_based_closed_to_internal(start: int, end: int) -> tuple[int, int]:
    """Convert a 1-based, closed interval to ``[start, end)``."""
    if start < 1 or end < start:
        raise CoordinateError("1-based closed coordinates require 1 <= start <= end")
    return start - 1, end


def internal_to_one_based_closed(start: int, end: int) -> tuple[int, int]:
    """Convert ``[start, end)`` to a 1-based, closed interval."""
    if start < 0 or end <= start:
        raise CoordinateError("internal coordinates require 0 <= start < end")
    return start + 1, end


def chromosome_alias(chrom: str) -> str:
    """Return a controlled canonical chromosome spelling.

    Only unambiguous aliases are changed; arbitrary fuzzy matching is not performed.
    """
    aliases = {"1": "chr1", "X": "chrX", "MT": "chrM", "M": "chrM"}
    return aliases.get(chrom, chrom)
