"""Reproducibility manifest generation."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

from .. import __version__
from ..model.figure_spec import FigureSpec
from .hashing import hash_file


def build_manifest(
    spec: FigureSpec, *, transformations: list[dict[str, Any]] | None = None
) -> dict[str, Any]:
    """Build deterministic provenance metadata for a figure specification."""
    spec_hash = hashlib.sha256(spec.to_json().encode()).hexdigest()
    inputs: list[dict[str, Any]] = []
    for track in spec.tracks:
        if track.source is not None and Path(track.source.path).is_file():
            source = Path(track.source.path)
            inputs.append(
                {"path": str(source), "size": source.stat().st_size, "sha256": hash_file(source)}
            )
    return {
        "software": {
            "name": "GenomeCanvas", "version": __version__, "python": platform.python_version()
        },
        "figure_spec_hash": spec_hash,
        "inputs": inputs,
        "region": spec.region.dict(),
        "transformations": transformations or [item.dict() for item in spec.transformations],
        "scale_decisions": [],
        "rendering": {"backend": "matplotlib"},
        "timestamp": None,
    }


def write_manifest(manifest: dict[str, Any], path: str | Path) -> None:
    """Write stable JSON provenance."""
    Path(path).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
