"""What the chat agent may do: the same actions a project editor has in the Studio UI.

Every tool runs as the signed-in user inside one project. Plots go through ``Studio.submit_job``,
so project access, per-user job limits, the queue, worker isolation and the audit log all apply.
The agent cannot read arbitrary files, run programs or see paths outside the project.
"""

from __future__ import annotations

import base64
import io
import json
import time
from pathlib import Path

from ..plan import PLOTS, job_body
from ..survey import survey

FIGURE_EDGE = 1400  # px; enough to judge a figure, small enough to keep image tokens modest

TOOLS = [
    dict(
        name="survey_data",
        description=(
            "Classify every file in this project (role, label, mark/factor, replicate, group, "
            "control), detect the genome, list problems, and return ranked figure "
            "recommendations with ready-to-run requests and the questions to ask the user. "
            "Call it first, before recommending anything."
        ),
        schema={"type": "object", "properties": {}, "additionalProperties": False},
    ),
    dict(
        name="read_skill",
        description=(
            "Read a GenomeCanvas guide: settings, defaults and pitfalls for one plot family, or "
            "the runner reference (request format and shared style settings). Read the family "
            "guide before changing settings beyond a drafted request."
        ),
        schema={
            "type": "object",
            "properties": {
                "name": {
                    "type": "string",
                    "enum": [
                        "enrichment",
                        "relationships",
                        "tracks",
                        "differential",
                        "transform",
                        "plot",
                    ],
                }
            },
            "required": ["name"],
            "additionalProperties": False,
        },
    ),
    dict(
        name="run_plot",
        description=(
            "Render one figure and show it to the user. `request` is a GenomeCanvas request: "
            '{"plot": "<slug>", "title": str, "files": [{"name" or "label", "role", "label", '
            '"group"}], "settings": {...}}. Refer to files by the names or labels that '
            "survey_data returned; never invent them. Returns status, warnings, the methods "
            "text and the rendered figure, which you should inspect before describing it. "
            f"Plot slugs: {', '.join(PLOTS)}."
        ),
        schema={
            "type": "object",
            "properties": {"request": {"type": "object"}},
            "required": ["request"],
            "additionalProperties": False,
        },
    ),
]


class Toolbox:
    def __init__(self, app, user, project, ip=None, emit=None):
        self.app, self.user, self.project, self.ip = app, user, project, ip
        self.emit = emit or (lambda event: None)
        self.skills = Path(app.config.static_dir).parent / "skills"

    def call(self, name, arguments):
        """Run a tool; returns (text, image_ref or None, is_error)."""
        try:
            if name == "survey_data":
                return self.survey_data(), None, False
            if name == "read_skill":
                return self.read_skill(arguments.get("name")), None, False
            if name == "run_plot":
                return self.run_plot(arguments.get("request"))
            return f"Unknown tool {name!r}", None, True
        except Exception as exc:  # noqa: BLE001 - every failure goes back to the model as text
            return str(getattr(exc, "args", [exc])[0] if exc.args else exc), None, True

    # ------------------------------------------------------------------ tools
    def files(self):
        return self.app.projects.files(self.project)

    def survey_data(self):
        files = self.files()
        names = {str(Path(f["path"]).resolve()): f for f in files}
        try:
            folder = self.app.projects.project_dir(self.project)
        except ValueError:
            folder = None
        # The project folder (with samples.tsv / project.yaml) plus uploads and shared files.
        paths = [
            f["path"] for f in files if not folder or not Path(f["path"]).is_relative_to(folder)
        ]
        report = survey(([str(folder)] if folder else []) + paths)

        def name_of(path):  # show project file names, never server paths
            row = names.get(str(Path(path).resolve()))
            return row["name"] if row else Path(path).name

        for item in report["files"]:
            item["path"] = name_of(item["path"])
        report["skipped"] = [dict(s, path=Path(s["path"]).name) for s in report["skipped"]]
        for rec in report["recommendations"]:
            rec.pop("request_file", None)
            for ref in (rec.get("request") or {}).get("files", []):
                ref["name"] = name_of(ref.pop("path"))
        for key in ("project", "sample_sheet"):
            if isinstance(report.get(key), dict):
                report[key].pop("path", None)
            elif report.get(key):
                report[key] = Path(report[key]).name
        return json.dumps(report, default=str)

    def read_skill(self, name):
        folder = "genomecanvas-plot" if name == "plot" else f"genomecanvas-{name}"
        path = self.skills / folder / "SKILL.md"
        if not path.is_file():
            raise ValueError(f"No guide named {name!r}")
        return path.read_text(encoding="utf-8")

    def run_plot(self, request):
        if not isinstance(request, dict):
            raise ValueError("run_plot needs a `request` object")
        body = job_body(request, self.files(), self.project.get("defaults") or {})
        job_id = self.app.submit_job(self.user, self.project, body, self.ip)
        self.emit(dict(type="job", job=job_id, plot=request.get("plot")))
        deadline = time.time() + self.app.config.job_timeout + 60
        while True:
            job = self.app.db.one("SELECT * FROM jobs WHERE id = ?", job_id)
            if job["status"] not in ("queued", "running"):
                break
            if time.time() > deadline:
                raise ValueError("The figure did not finish in time")
            time.sleep(0.5)
        if job["status"] != "complete":
            return f"The figure failed: {job['error']}", None, True
        folder = self.app.jobs.job_dir(self.project["id"], job_id)
        self.emit(dict(type="figure", job=job_id, url=f"/results/{job_id}/figure.png"))
        methods = folder / "methods.txt"
        result = job["result"] or {}
        text = json.dumps(
            dict(
                status="complete",
                job=job_id,
                elapsed=round(job["elapsed"] or 0, 1),
                warnings=result.get("warnings", []),
                methods=methods.read_text(encoding="utf-8").strip() if methods.is_file() else "",
                note="The user can already see this figure in the chat.",
            )
        )
        return text, dict(job=job_id, project=self.project["id"]), False


def load_image(app, ref):
    """Downscaled PNG of a job's figure as an image block payload, or None if it has expired."""
    from PIL import Image

    path = app.jobs.job_dir(ref["project"], ref["job"]) / "figure.png"
    if not path.is_file():
        return None
    with Image.open(path) as image:
        image = image.convert("RGB")
        image.thumbnail((FIGURE_EDGE, FIGURE_EDGE))
        buffer = io.BytesIO()
        image.save(buffer, format="PNG", optimize=True)
    return dict(media_type="image/png", data=base64.b64encode(buffer.getvalue()).decode())
