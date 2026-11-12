"""Plot requests as agents and scripts write them: stable slugs, UI defaults and file matching.

Shared by ``skills/genomecanvas-plot/scripts/gcplot.py`` and the Studio's chat agent, so a request
renders the same figure from the command line, a script or the chat window.
"""

from __future__ import annotations

from pathlib import Path

# Stable slugs for CATALOG[category][1][plot] in genomecanvas.studio.engine.
PLOTS = {
    "tracks/stacked": (0, 0),
    "tracks/overlay": (0, 1),
    "tracks/genes": (0, 2),
    "enrichment/heatmap": (1, 0),
    "enrichment/profile": (1, 1),
    "enrichment/profile-heatmap": (1, 2),
    "enrichment/clustered-heatmap": (1, 3),
    "relationships/correlation": (2, 0),
    "relationships/pca": (2, 1),
    "relationships/pairgrid": (2, 2),
    "relationships/sample-qc": (2, 3),
    "differential/volcano": (3, 0),
    "differential/peak-center": (3, 1),
    "transform/difference": (4, 0),
    "transform/ratio": (4, 1),
    "transform/log2-ratio": (4, 2),
    "transform/background": (4, 3),
}
SKILLS = ["tracks", "enrichment", "relationships", "differential", "transform"]


def ui_defaults(category, plot, settings):
    """Mirror effectiveSettings() in web/app.js where the UI default differs from the engine's."""
    out = dict(settings)
    out.setdefault("size_unit", "mm")
    if category == 1 or (category == 3 and plot == 1):
        out.setdefault("flank", 1000)
        out.setdefault("binsize", 50)
    if category == 1:
        out.setdefault("fill_alpha", 0.18)
    if category == 4 and plot == 3:
        # The UI sends its "Background flank" field as `flank` for this plot only.
        out["flank"] = out.get("bg_flank") or out.get("flank") or 5000
    if category == 2:
        out.setdefault("rel_stat", "sum" if out.get("rel_mode", "regions") == "regions" else "mean")
    return out


def resolve_plot(request):
    """Return (category, plot) from a `plot` slug or explicit integers."""
    slug = request.get("plot")
    if isinstance(slug, str):
        if slug not in PLOTS:
            raise ValueError(f"Unknown plot {slug!r}; choose one of {', '.join(PLOTS)}")
        return PLOTS[slug]
    if "category" in request:
        return int(request["category"]), int(request.get("plot") or 0)
    raise ValueError("The request needs `plot` (e.g. 'relationships/pca') or `category` + `plot`")


def slug_of(category, plot):
    return next((s for s, cp in PLOTS.items() if cp == (category, plot)), f"{category}/{plot}")


def match_file(ref, files):
    """Project file id for a reference by id, file name, path basename or label (in that order)."""
    if isinstance(ref, dict) and ref.get("id"):
        return ref["id"]
    key = ref.get("name") or ref.get("path") or ref.get("label") if isinstance(ref, dict) else ref
    key = str(key or "")
    for field in ("id", "name", "label"):
        hits = [
            f for f in files if f[field] == key or (field == "name" and f[field] == Path(key).name)
        ]
        if len(hits) == 1:
            return hits[0]["id"]
        if len(hits) > 1:
            raise ValueError(f"{key!r} matches several project files; use the file id")
    raise ValueError(f"{key!r} is not a file in this project")


def job_body(request, files, project_defaults=None):
    """The Studio jobs-API body for a request; `files` are the project's files (id/name/label)."""
    category, plot = resolve_plot(request)
    # Project defaults sit under the request, as when the UI opens the project.
    settings = ui_defaults(
        category, plot, {**(project_defaults or {}), **(request.get("settings") or {})}
    )
    for key in ("rel_blacklist", "qc_sheet", "qc_matrix"):
        if settings.get(key):
            settings[key] = match_file({"name": settings[key]}, files)
    items = []
    for ref in request.get("files") or []:
        item = {k: v for k, v in (ref if isinstance(ref, dict) else {}).items() if k != "path"}
        item["id"] = match_file(ref, files)
        items.append(item)
    return dict(
        category=category,
        plot=plot,
        files=items,
        settings=settings,
        title=request.get("title") or settings.get("title"),
        palette=settings.get("palette"),
        advanced=request.get("advanced") or {},
    )
