"""Deterministic greedy packing of genomic intervals into display rows."""

from __future__ import annotations

from collections.abc import Iterable

from ..core.errors import CoordinateError


def pack_intervals(intervals: Iterable[tuple[int, int]]) -> tuple[tuple[int, ...], ...]:
    """Assign intervals to the first compatible row in deterministic order.

    Intervals are sorted by ``(start, end, original_position)``. Touching
    intervals are compatible because genomic intervals are half-open.
    Returned rows contain original input indices, making the result useful for
    rendering without changing feature order or coordinates.
    """
    indexed = list(enumerate(intervals))
    for index, (start, end) in indexed:
        if start < 0 or end <= start:
            raise CoordinateError(f"interval {index} must satisfy 0 <= start < end")
    ordered = sorted(indexed, key=lambda item: (item[1][0], item[1][1], item[0]))
    row_ends: list[int] = []
    rows: list[list[int]] = []
    for index, (_, end) in ordered:
        start = indexed[index][1][0]
        for row, row_end in enumerate(row_ends):
            if start >= row_end:
                rows[row].append(index)
                row_ends[row] = end
                break
        else:
            rows.append([index])
            row_ends.append(end)
    return tuple(tuple(row) for row in rows)
