"""One entry point for every GenomeCanvas Studio plot, for local use and for agents.

    python3 skills/genomecanvas-plot/scripts/gcplot.py list
    python3 skills/genomecanvas-plot/scripts/gcplot.py survey DIR_OR_FILES ... [--write DIR]
    python3 skills/genomecanvas-plot/scripts/gcplot.py inspect FILE [FILE ...]
    python3 skills/genomecanvas-plot/scripts/gcplot.py run REQUEST.json [--out DIR] [--set k=v ...]
    python3 skills/genomecanvas-plot/scripts/gcplot.py run REQUEST.json --server URL --project ID

`run` without --server calls the Studio engine in-process on local file paths. With --server it
submits the same request to a running Studio (files referenced by id, file name or label) and
downloads the results. Both modes apply the Studio UI defaults, so a request renders the same figure
either way. Every command prints one JSON object on stdout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from genomecanvas.studio.plan import PLOTS, SKILLS, job_body, slug_of, ui_defaults  # noqa: E402
from genomecanvas.studio.plan import resolve_plot as _resolve_plot  # noqa: E402


def fail(message, code=1):
    print(json.dumps({"status": "error", "error": str(message)}, ensure_ascii=False))
    sys.exit(code)


def emit(obj):
    print(json.dumps(obj, indent=2, ensure_ascii=False, default=str))


def resolve_plot(request):
    """Return (category, plot) from a `plot` slug or explicit integers."""
    try:
        return _resolve_plot(request)
    except ValueError as exc:
        fail(str(exc))

def parse_sets(pairs):
    settings = {}
    for pair in pairs or []:
        key, _, raw = pair.partition("=")
        if not key or not _:
            fail(f"--set expects key=value, got {pair!r}")
        try:
            settings[key] = json.loads(raw)
        except json.JSONDecodeError:
            settings[key] = raw
    return settings


def load_request(path, sets):
    try:
        request = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"Cannot read request {path}: {exc}")
    request.setdefault("settings", {})
    request["settings"].update(parse_sets(sets))
    request.setdefault("files", [])
    # Sample QC reads its sheet and matrix from settings, not from `files`.
    if not request["files"] and request.get("plot") != "relationships/sample-qc":
        fail("The request needs a non-empty `files` list")
    return request


# ---------------------------------------------------------------- local mode


def local_record(item, base):
    from genomecanvas.studio.engine import SUFFIXES, inspect_file

    path = Path(item["path"]).expanduser()
    if not path.is_absolute():
        path = (base / path).resolve()
    if not path.is_file():
        fail(f"Input file not found: {path}")
    suffixes = [s.lower() for s in path.suffixes]
    guess = next((SUFFIXES[s] for s in reversed(suffixes) if s in SUFFIXES), None)
    role = item.get("role") or guess
    if role not in ("signal", "regions", "genes", "table"):
        fail(f"Set `role` (signal, regions, genes or table) for {path.name}")
    return dict(
        id=hashlib.sha1(str(path).encode()).hexdigest()[:32],
        path=str(path),
        name=path.name,
        role=role,
        source="upload",
        group=str(item.get("group") or "group"),
        label=str(item.get("label") or path.name.split(".")[0]),
        color=str(item.get("color") or ""),
        cmap=str(item.get("cmap") or ""),
        track=str(item.get("track") or ""),
        info=inspect_file(path, role),
    )


def run_local(request, out, base, cache_dir):
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir.parent / "mpl"))
    from genomecanvas.studio.engine import run_plot

    category, plot = resolve_plot(request)
    settings = ui_defaults(category, plot, request.get("settings") or {})
    records = [local_record(item, base) for item in request["files"]]
    for key in ("qc_sheet", "qc_matrix"):
        if (category, plot) == (2, 3) and settings.get(key):
            records.append(
                {**local_record({"path": settings[key], "role": "table"}, base), "role": key}
            )
    blacklist = settings.get("rel_blacklist")
    if blacklist:
        # Same as the server: authorise the file, then tag it so it filters bins but is not plotted.
        records.append(
            {**local_record({"path": blacklist, "role": "regions"}, base), "role": "blacklist"}
        )
    payload = dict(
        category=category,
        plot=plot,
        settings=settings,
        title=request.get("title") or settings.get("title"),
        palette=settings.get("palette"),
        advanced=request.get("advanced") or {},
    )
    out.mkdir(parents=True, exist_ok=True)
    began = time.time()
    try:
        warnings, counts = run_plot(payload, records, out, cache_dir)
    except Exception as exc:  # noqa: BLE001 - reported as JSON for the caller
        fail(exc)
    return dict(warnings=warnings, counts=counts, elapsed=round(time.time() - began, 2))


# ---------------------------------------------------------------- server mode


class Studio:
    def __init__(self, url, cookie):
        self.url = url.rstrip("/")
        self.headers = {"X-GenomeCanvas": "studio", "Origin": self.url}
        if cookie:
            self.headers["Cookie"] = cookie if "=" in cookie else f"gc_session={cookie}"

    def call(self, path, body=None):
        data = json.dumps(body).encode() if body is not None else None
        headers = {**self.headers, **({"Content-Type": "application/json"} if data else {})}
        request = urllib.request.Request(self.url + path, data=data, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.load(response)
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode(errors="replace")
            try:
                detail = json.loads(detail).get("error", detail)
            except (json.JSONDecodeError, AttributeError):
                pass
            fail(f"{exc.code} on {path}: {detail}")
        except urllib.error.URLError as exc:
            fail(f"Cannot reach {self.url}: {exc.reason}")

    def download(self, path, target):
        request = urllib.request.Request(self.url + path, headers=self.headers)
        with urllib.request.urlopen(request, timeout=120) as response:
            target.write_bytes(response.read())


def run_server(request, out, studio, project_id, poll):
    detail = studio.call(f"/api/projects/{urllib.parse.quote(project_id)}")
    try:
        body = job_body(request, detail["files"], detail["project"].get("defaults", {}))
    except ValueError as exc:
        fail(str(exc))
    job = studio.call(f"/api/projects/{urllib.parse.quote(project_id)}/jobs", body)
    deadline = time.time() + 900
    while True:
        state = studio.call(f"/api/jobs/{job['id']}")
        if state["status"] == "complete":
            break
        if state["status"] not in ("queued", "running"):
            fail(state.get("error") or f"Job ended with status {state['status']}")
        if time.time() > deadline:
            fail(f"Job {job['id']} did not finish within 15 minutes")
        time.sleep(poll)
    out.mkdir(parents=True, exist_ok=True)
    for item in state.get("files", []):
        studio.download(item["url"], out / item["name"])
    return dict(
        job=job["id"],
        warnings=state.get("warnings", []),
        counts=state.get("counts"),
        elapsed=state.get("elapsed"),
        result_urls=[studio.url + f["url"] for f in state.get("files", [])],
    )


# ---------------------------------------------------------------- commands


def cmd_list(_args):
    from genomecanvas.studio.engine import CATALOG, REQUIREMENTS

    plots = []
    for slug, (c, p) in PLOTS.items():
        plots.append(
            dict(
                plot=slug,
                name=CATALOG[c][1][p],
                category=CATALOG[c][0],
                needs=REQUIREMENTS[c][p],
                skill=f"skills/genomecanvas-{SKILLS[c]}/SKILL.md",
            )
        )
    emit(dict(status="ok", plots=plots))


def cmd_inspect(args):
    from genomecanvas.studio.engine import SUFFIXES, inspect_file

    rows = []
    for raw in args.files:
        path = Path(raw).expanduser().resolve()
        suffixes = [s.lower() for s in path.suffixes]
        role = next((SUFFIXES[s] for s in reversed(suffixes) if s in SUFFIXES), None)
        info = inspect_file(path, role or "regions") if path.is_file() else {"error": "not found"}
        rows.append(dict(path=str(path), role=role, info=info))
    emit(dict(status="ok", files=rows))


def cmd_survey(args):
    from genomecanvas.studio.survey import survey

    emit(survey(args.paths, write_dir=args.write))


def cmd_files(args):
    detail = Studio(args.server, args.cookie).call(
        f"/api/projects/{urllib.parse.quote(args.project)}"
    )
    emit(
        dict(
            status="ok",
            project=detail["project"].get("name"),
            defaults=detail["project"].get("defaults", {}),
            files=[
                {k: f[k] for k in ("id", "name", "role", "label", "group")} for f in detail["files"]
            ],
        )
    )


def cmd_run(args):
    request = load_request(args.request, args.set)
    category, plot = resolve_plot(request)
    slug = slug_of(category, plot).replace("/", "-")
    out = Path(args.out or Path.cwd() / "gcanvas_out" / f"{slug}-{time.strftime('%Y%m%d-%H%M%S')}")
    out = out.expanduser().resolve()
    if args.server:
        result = run_server(request, out, Studio(args.server, args.cookie), args.project, args.poll)
    else:
        cache = Path(args.cache or Path.home() / ".cache" / "genomecanvas" / "skills").expanduser()
        result = run_local(request, out, Path(args.request).resolve().parent, cache)
    methods = out / "methods.txt"
    emit(
        dict(
            status="complete",
            plot=plot_name(category, plot),
            out=str(out),
            figure=str(out / "figure.png"),
            files=sorted(p.name for p in out.iterdir() if p.is_file()),
            methods=methods.read_text(encoding="utf-8").strip() if methods.is_file() else None,
            **result,
        )
    )


def plot_name(category, plot):
    from genomecanvas.studio.engine import CATALOG

    return CATALOG[category][1][plot]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="plots, their slugs and minimum inputs")
    survey = sub.add_parser(
        "survey", help="classify a data folder, flag problems and recommend figures"
    )
    survey.add_argument("paths", nargs="+", help="folders (searched recursively) or files")
    survey.add_argument("--write", help="also write each recommended request as JSON here")
    inspect = sub.add_parser("inspect", help="role, chromosomes and extent of input files")
    inspect.add_argument("files", nargs="+")
    run = sub.add_parser("run", help="render one figure from a request JSON")
    run.add_argument("request")
    run.add_argument("--out", help="output directory (default ./gcanvas_out/<plot>-<time>)")
    run.add_argument("--set", action="append", metavar="KEY=VALUE", help="override one setting")
    run.add_argument("--cache", help="matrix/summary cache for local runs")
    for p in (run, sub.add_parser("files", help="list a Studio project's files (server mode)")):
        p.add_argument("--server", help="Studio base URL, e.g. http://127.0.0.1:8765")
        p.add_argument("--project", default="local", help="Studio project id (server mode)")
        p.add_argument(
            "--cookie",
            default=os.environ.get("GCANVAS_SESSION"),
            help="session cookie for sign-in deployments (or env GCANVAS_SESSION)",
        )
    run.add_argument("--poll", type=float, default=0.5, help="seconds between job status checks")
    args = parser.parse_args()
    if args.command == "files" and not args.server:
        fail("files needs --server")
    commands = dict(
        list=cmd_list, survey=cmd_survey, inspect=cmd_inspect, run=cmd_run, files=cmd_files
    )
    commands[args.command](args)


if __name__ == "__main__":
    main()
