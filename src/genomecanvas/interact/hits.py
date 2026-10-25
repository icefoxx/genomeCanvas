"""Hit-target registry recorded while drawing, for later pointer mapping."""

from __future__ import annotations

from ..model.overlays import HitTarget


class HitRegistry:
    """Collect pickable figure objects during a render pass."""

    def __init__(self) -> None:
        self.targets: list[HitTarget] = []

    def add(self, target: HitTarget) -> HitTarget:
        """Record one hit target and return it."""
        self.targets.append(target)
        return target

    def to_json(self) -> list[dict[str, object]]:
        """Serialize targets for a sidecar ``.hits.json`` file."""
        return [item.dict() for item in self.targets]


def hit_at_genome(
    targets: list[HitTarget],
    chrom: str,
    position: int,
    *,
    track_id: str | None = None,
    axes_index: int | None = None,
) -> HitTarget | None:
    """Return the last target covering ``position`` on ``chrom``.

    Interactive backends can map a click in data coordinates to this query.
    Later targets win so overlays sit above tracks.
    """
    match: HitTarget | None = None
    for target in targets:
        if target.chrom != chrom:
            continue
        if not (target.start <= position < target.end or target.start == target.end == position):
            if not (target.start <= position <= max(target.end, target.start)):
                continue
        if track_id is not None and target.track_id != track_id:
            continue
        if axes_index is not None and target.axes_index != axes_index:
            continue
        match = target
    return match
