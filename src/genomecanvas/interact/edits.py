"""Apply FigureEdit payloads to a FigureSpec without re-parsing YAML."""

from __future__ import annotations

from collections.abc import Callable

from ..model.figure_spec import FigureSpec, TrackSpec
from ..model.overlays import FigureEdit, OverlayAnnotation
from ..model.style import FigureStyle


def overlay_by_id(spec: FigureSpec, overlay_id: str) -> OverlayAnnotation | None:
    """Return the overlay with ``overlay_id``, if present."""
    for item in spec.overlays:
        if item.id == overlay_id:
            return item
    return None


def apply_edit(spec: FigureSpec, edit: FigureEdit) -> FigureSpec:
    """Return a new FigureSpec with ``edit`` applied.

    Click handlers should convert pointer events into ``FigureEdit`` values and
    call this function. Rendering stays a pure function of the spec.
    """
    if edit.action == "add_overlay":
        if edit.overlay is None:
            raise ValueError("add_overlay requires edit.overlay")
        if overlay_by_id(spec, edit.overlay.id):
            raise ValueError(f"overlay {edit.overlay.id!r} already exists")
        return spec.copy(update={"overlays": [*spec.overlays, edit.overlay]})
    if edit.action == "update_overlay":
        if edit.overlay is None:
            raise ValueError("update_overlay requires edit.overlay")
        overlays = [
            edit.overlay if item.id == edit.overlay.id else item for item in spec.overlays
        ]
        if all(item.id != edit.overlay.id for item in spec.overlays):
            overlays.append(edit.overlay)
        return spec.copy(update={"overlays": overlays})
    if edit.action == "remove_overlay":
        target = edit.overlay_id or (edit.overlay.id if edit.overlay else None)
        if not target:
            raise ValueError("remove_overlay requires overlay_id")
        return spec.copy(
            update={"overlays": [item for item in spec.overlays if item.id != target]}
        )
    if edit.action == "update_track_style":
        return _map_track(
            spec, edit.track_id, lambda track: _merge_style(track, edit.style or {})
        )
    if edit.action == "set_track_height":
        if edit.height is None or edit.height <= 0:
            raise ValueError("set_track_height requires a positive height")
        return _map_track(
            spec, edit.track_id, lambda track: track.copy(update={"height": edit.height})
        )
    if edit.action == "set_track_visibility":
        if edit.visibility is None:
            raise ValueError("set_track_visibility requires visibility")
        return _map_track(
            spec,
            edit.track_id,
            lambda track: track.copy(update={"visibility": edit.visibility}),
        )
    if edit.action == "update_figure_style":
        current = spec.style.dict()
        current.update(edit.style or {})
        return spec.copy(update={"style": FigureStyle.parse_obj(current)})
    if edit.action == "set_title":
        style = spec.style.copy(update={"title": edit.title})
        return spec.copy(update={"style": style})
    raise ValueError(f"unsupported edit action {edit.action!r}")


def _map_track(
    spec: FigureSpec, track_id: str | None, mutate: Callable[[TrackSpec], TrackSpec]
) -> FigureSpec:
    if not track_id:
        raise ValueError("track_id is required for this edit")
    tracks = []
    found = False
    for track in spec.tracks:
        if track.id == track_id:
            tracks.append(mutate(track))
            found = True
        else:
            tracks.append(track)
    if not found:
        raise ValueError(f"unknown track {track_id!r}")
    return spec.copy(update={"tracks": tracks})


def _merge_style(track: TrackSpec, style: dict[str, object]) -> TrackSpec:
    merged = dict(track.style)
    merged.update(style)
    return track.copy(update={"style": merged})
