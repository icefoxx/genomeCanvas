"""Projects: the isolation unit. Deliverables come from a manifest; uploads are per project."""

from __future__ import annotations

import csv
import re
import shutil
import time
from pathlib import Path

import yaml

from .db import Database, new_id
from .engine import SUFFIXES, inspect_file

SHARED = "_shared"
MANIFESTS = ("samples.tsv", "genomecanvas.tsv", "manifest.tsv")
PROJECT_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{1,63}$")
ROLE_ALIASES = {
    "signal": "signal",
    "bigwig": "signal",
    "bw": "signal",
    "bedgraph": "signal",
    "regions": "regions",
    "region": "regions",
    "bed": "regions",
    "peaks": "regions",
    "peak": "regions",
    "genes": "genes",
    "gene": "genes",
    "gtf": "genes",
    "annotation": "genes",
}
FIXTURES = Path(__file__).resolve().parents[3] / "tests" / "fixtures"
DEMOS = {
    "demo": dict(
        name="Demo · H1 ESC OCT4 locus",
        description="Public trackplot demo: 7 ChIP-seq BigWigs on hg19 chr6 with peaks, chromHMM and refGene.",
        genome="hg19",
        defaults=dict(
            chrom="chr6",
            start=31125000,
            end=31145000,
            assembly="hg19",
            flank=1000,
            binsize=50,
            alignment="center",
            xlab="distance from region center",
            colorbar_label="signal",
            numerator="TF",
            denominator="Histone",
            title="H1 ESC · OCT4 locus",
        ),
        files=[
            ("trackplot_demo/bw/H1_k4me3.bw", "signal", "Histone", "H3K4me3"),
            ("trackplot_demo/bw/H1_k27ac.bw", "signal", "Histone", "H3K27ac"),
            ("trackplot_demo/bw/H1_k4me1.bw", "signal", "Histone", "H3K4me1"),
            ("trackplot_demo/bw/H1_H2az.bw", "signal", "Histone", "H2A.Z"),
            ("trackplot_demo/bw/H1_Ctcf.bw", "signal", "TF", "CTCF"),
            ("trackplot_demo/bw/H1_Oct4.bw", "signal", "TF", "OCT4"),
            ("trackplot_demo/bw/H1_Nanog.bw", "signal", "TF", "NANOG"),
            ("trackplot_demo/narrowpeak/H1_chromHMM.bed", "regions", "chromHMM", "chromHMM states"),
            ("trackplot_demo/narrowpeak/H1_k4me1.bed", "regions", "H3K4me1 peaks", "H3K4me1 peaks"),
            ("trackplot_demo/chr6_refGene.gtf", "genes", "annotation", "refGene"),
        ],
    ),
    "demo-ce10": dict(
        name="Demo · SeqPlots C. elegans",
        description="SeqPlots inst/extdata: one H3K4me3 BigWig, peak calls and transcripts (ce10 chrI 0–100 kb).",
        genome="ce10",
        defaults=dict(
            chrom="chrI",
            start=0,
            end=100000,
            assembly="ce10",
            flank=500,
            binsize=25,
            alignment="center",
            xlab="distance from center",
            colorbar_label="coverage",
            title="SeqPlots example",
        ),
        files=[
            ("seqplots_demo/GSM1208360_chrI_100Kb_q5_sample.bw", "signal", "H3K4me3", "H3K4me3"),
            ("seqplots_demo/GSM1208361_chrI_100Kb_PeakCalls.bed", "regions", "Peaks", "Peaks"),
            (
                "seqplots_demo/Transcripts_ce10_chrI_100Kb.bed",
                "regions",
                "Transcripts",
                "Transcripts",
            ),
            ("seqplots_demo/transcripts_ce10_chrI_100Kb.gtf", "genes", "annotation", "Genes"),
        ],
    ),
}


class Projects:
    def __init__(self, db: Database, config):
        self.db = db
        self.config = config
        self.data_root = Path(config.data_root).expanduser().resolve() if config.data_root else None
        self.shared_root = (
            Path(config.shared_root).expanduser().resolve() if config.shared_root else None
        )
        self.uploads = config.workspace_path / "uploads"

    # ------------------------------------------------------------ access
    def is_demo(self, project_id):
        return self.config.demo_project and project_id in DEMOS

    def role(self, user, project_id):
        project = self.db.one("SELECT * FROM projects WHERE id = ?", project_id)
        if not project or user is None:
            return None, project
        role = self.db.project_role(user, project_id)
        if role is None and self.is_demo(project_id):
            role = "editor"
        if role and role != "admin" and project["status"] != "active":
            return None, project
        if (
            role
            and role != "admin"
            and project.get("expires")
            and project["expires"] < time.strftime("%Y-%m-%d")
        ):
            return None, project
        return role, project

    def accessible(self, user):
        projects = {p["id"]: p for p in self.db.projects_for(user)}
        if self.config.demo_project:
            for pid in DEMOS:
                if pid not in projects:
                    project = self.db.one("SELECT * FROM projects WHERE id = ?", pid)
                    if project:
                        projects[pid] = {**project, "role": "editor"}
        today = time.strftime("%Y-%m-%d")
        out = []
        for project in projects.values():
            if project["role"] != "admin" and (
                project["status"] != "active"
                or (project.get("expires") and project["expires"] < today)
            ):
                continue
            project["demo"] = project["id"] in DEMOS
            out.append(project)
        return sorted(out, key=lambda p: (p["demo"], -(p["created"] or 0)))

    # ------------------------------------------------------------ admin operations
    def create(
        self,
        project_id,
        name=None,
        data_dir="",
        genome="",
        description="",
        expires=None,
        quota_gb=None,
    ):
        if not PROJECT_ID.match(project_id or ""):
            raise ValueError("Project IDs use letters, digits, '.', '_' or '-' (2–64 characters)")
        if self.db.one("SELECT id FROM projects WHERE id = ?", project_id):
            raise ValueError(f"Project {project_id} already exists")
        if data_dir:
            self.project_dir(dict(data_dir=data_dir))
        quota = quota_gb if quota_gb is not None else self.config.default_quota_gb
        self.db.insert(
            "projects",
            id=project_id,
            name=name or project_id,
            description=description,
            data_dir=data_dir,
            genome=genome,
            status="active",
            expires=expires or None,
            quota_bytes=int(quota * 1024**3),
            defaults={},
            created=time.time(),
        )
        return self.db.one("SELECT * FROM projects WHERE id = ?", project_id)

    def project_dir(self, project) -> Path | None:
        rel = project.get("data_dir") or ""
        if not rel:
            return None
        if self.data_root is None:
            raise ValueError(
                "Set data_root in the studio configuration before attaching project data"
            )
        path = (self.data_root / rel).resolve()
        if not path.is_relative_to(self.data_root):
            raise ValueError("Project data must live under data_root")
        if not path.is_dir():
            raise ValueError(f"Project data folder not found: {path}")
        return path

    def set_member(self, project_id, user_id, role):
        if role not in {"editor", "viewer", None}:
            raise ValueError("Role must be editor or viewer")
        self.db.run("DELETE FROM members WHERE project_id = ? AND user_id = ?", project_id, user_id)
        if role:
            self.db.insert("members", project_id=project_id, user_id=user_id, role=role)

    # ------------------------------------------------------------ scanning
    def _upsert(
        self,
        project_id,
        path,
        *,
        role,
        label,
        group,
        color="",
        cmap="",
        track="",
        source,
        sort=0,
        existing,
    ):
        path = Path(path)
        stat = path.stat()
        old = existing.pop(str(path), None)
        info = (
            old["info"]
            if old and old["size"] == stat.st_size and old["mtime"] == stat.st_mtime
            else None
        )
        if info is None or (old and old["role"] != role):
            info = inspect_file(path, role)
        values = dict(
            role=role,
            label=label,
            grp=group,
            color=color,
            cmap=cmap,
            track=track,
            source=source,
            size=stat.st_size,
            mtime=stat.st_mtime,
            sort=sort,
            info=info,
            name=path.name,
        )
        if old:
            self.db.update("files", "id", old["id"], **values)
            return old["id"]
        file_id = new_id()
        self.db.insert(
            "files",
            id=file_id,
            project_id=project_id,
            path=str(path),
            created=time.time(),
            **values,
        )
        return file_id

    def scan(self, project_id):
        """Synchronise a project's deliverables with its folder and manifest."""
        project = self.db.one("SELECT * FROM projects WHERE id = ?", project_id)
        if not project:
            raise ValueError(f"No project {project_id}")
        if project_id in DEMOS:
            return self._scan_demo(project_id)
        folder = self.project_dir(project)
        existing = {
            f["path"]: f
            for f in self.db.all(
                "SELECT * FROM files WHERE project_id = ? AND source IN ('manifest', 'scan')",
                project_id,
            )
        }
        report = dict(added=0, updated=0, removed=0, warnings=[])
        if folder is None:
            return report
        settings = self._project_yaml(project, folder, report)
        manifest = next((folder / m for m in MANIFESTS if (folder / m).is_file()), None)
        rows = (
            self._manifest_rows(folder, manifest, report) if manifest else self._folder_rows(folder)
        )
        known = set(existing)
        for order, row in enumerate(rows):
            self._upsert(project_id, sort=order, existing=existing, **row)
        for stale in existing.values():
            self.db.run("DELETE FROM files WHERE id = ?", stale["id"])
        added = sum(str(row["path"]) not in known for row in rows)
        report.update(
            files=len(rows), added=added, updated=len(rows) - added, removed=len(existing)
        )
        if settings:
            self.db.update("projects", "id", project_id, **settings)
        return report

    def _project_yaml(self, project, folder, report):
        path = folder / "project.yaml"
        if not path.is_file():
            return {}
        try:
            data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except yaml.YAMLError as exc:
            report["warnings"].append(f"project.yaml: {exc}")
            return {}
        updates = {}
        for key in ("name", "description", "genome"):
            if data.get(key):
                updates[key] = str(data[key])
        if isinstance(data.get("defaults"), dict):
            updates["defaults"] = data["defaults"]
        self.db.run(
            "DELETE FROM templates WHERE project_id = ? AND created_by = 'manifest'", project["id"]
        )
        for item in data.get("templates") or []:
            if isinstance(item, dict) and item.get("name") and isinstance(item.get("config"), dict):
                self.db.insert(
                    "templates",
                    id=new_id(),
                    project_id=project["id"],
                    name=str(item["name"]),
                    description=str(item.get("description") or ""),
                    config=item["config"],
                    created_by="manifest",
                    created=time.time(),
                )
        return updates

    def _manifest_rows(self, folder, manifest, report):
        rows = []
        with manifest.open(encoding="utf-8") as handle:
            lines = [line for line in handle if line.strip() and not line.startswith("#")]
        for number, row in enumerate(csv.DictReader(lines, delimiter="\t"), 2):
            row = {str(k).strip().lower(): (v or "").strip() for k, v in row.items() if k}
            rel = row.get("file") or row.get("path")
            if not rel:
                report["warnings"].append(f"{manifest.name} line {number}: missing 'file' column")
                continue
            path = (folder / rel).resolve()
            if not path.is_relative_to(folder) or not path.is_file():
                report["warnings"].append(
                    f"{manifest.name} line {number}: {rel} not found in project folder"
                )
                continue
            role = ROLE_ALIASES.get(row.get("role", "").lower()) or SUFFIXES.get(
                path.suffix.lower()
            )
            if not role:
                report["warnings"].append(
                    f"{manifest.name} line {number}: unsupported file type {rel}"
                )
                continue
            rows.append(
                dict(
                    path=path,
                    role=role,
                    label=row.get("label") or path.name.split(".")[0],
                    group=row.get("group") or row.get("condition") or path.parent.name,
                    color=row.get("color", ""),
                    cmap=row.get("cmap", ""),
                    track=row.get("track", "").lower(),
                    source="manifest",
                )
            )
        return rows

    def _folder_rows(self, folder):
        rows = []
        for path in sorted(folder.rglob("*")):
            role = SUFFIXES.get(path.suffix.lower())
            if (
                not role
                or not path.is_file()
                or any(p.startswith(".") for p in path.relative_to(folder).parts)
            ):
                continue
            rows.append(
                dict(
                    path=path,
                    role=role,
                    label=path.name.split(".")[0],
                    group=path.parent.name,
                    source="scan",
                )
            )
            if len(rows) >= 500:
                break
        return rows

    def _scan_demo(self, project_id):
        spec = DEMOS[project_id]
        existing = {
            f["path"]: f
            for f in self.db.all("SELECT * FROM files WHERE project_id = ?", project_id)
        }
        if project_id == "demo-ce10":
            _transcripts_gtf(
                FIXTURES / "seqplots_demo/Transcripts_ce10_chrI_100Kb.bed",
                FIXTURES / "seqplots_demo/transcripts_ce10_chrI_100Kb.gtf",
            )
        count = 0
        for order, (rel, role, group, label) in enumerate(spec["files"]):
            path = FIXTURES / rel
            if path.is_file():
                self._upsert(
                    project_id,
                    path,
                    role=role,
                    label=label,
                    group=group,
                    source="example",
                    sort=order,
                    existing=existing,
                )
                count += 1
        return dict(files=count, warnings=[] if count else ["demo fixtures are missing"])

    def ensure_demos(self):
        if not self.config.demo_project:
            return
        for pid, spec in DEMOS.items():
            if not self.db.one("SELECT id FROM projects WHERE id = ?", pid):
                self.db.insert(
                    "projects",
                    id=pid,
                    name=spec["name"],
                    description=spec["description"],
                    data_dir="",
                    genome=spec["genome"],
                    status="active",
                    expires=None,
                    quota_bytes=int(1 * 1024**3),
                    defaults=spec["defaults"],
                    created=0,
                )
            self._scan_demo(pid)

    def scan_shared(self):
        """Register <shared_root>/<genome>/* annotation files, visible to every project."""
        if not self.shared_root or not self.shared_root.is_dir():
            return dict(files=0)
        existing = {
            f["path"]: f for f in self.db.all("SELECT * FROM files WHERE project_id = ?", SHARED)
        }
        count = 0
        for genome_dir in sorted(p for p in self.shared_root.iterdir() if p.is_dir()):
            for path in sorted(genome_dir.iterdir()):
                role = SUFFIXES.get(path.suffix.lower())
                if role and path.is_file():
                    self._upsert(
                        SHARED,
                        path,
                        role=role,
                        label=f"{genome_dir.name} · {path.stem}",
                        group=genome_dir.name,
                        source="shared",
                        sort=count,
                        existing=existing,
                    )
                    count += 1
        for stale in existing.values():
            self.db.run("DELETE FROM files WHERE id = ?", stale["id"])
        return dict(files=count)

    # ------------------------------------------------------------ files
    def files(self, project):
        rows = self.db.all(
            "SELECT * FROM files WHERE project_id = ? ORDER BY source = 'upload', sort, created",
            project["id"],
        )
        genome = (project.get("genome") or "").lower()
        shared = self.db.all("SELECT * FROM files WHERE project_id = ? ORDER BY sort", SHARED)
        rows += [f for f in shared if not genome or (f["grp"] or "").lower() == genome]
        return rows

    def resolve(self, project, items):
        """Authorise requested file ids for this project and merge per-request overrides."""
        if not isinstance(items, list) or not items:
            raise ValueError("Select at least one input file")
        if len(items) > 48:
            raise ValueError("Choose at most 48 input files")
        allowed = {f["id"]: f for f in self.files(project)}
        records = []
        for item in items:
            row = allowed.get(item.get("id") if isinstance(item, dict) else None)
            if row is None:
                raise ValueError(
                    "A selected file is not available in this project; refresh the page"
                )
            role = item.get("role") or row["role"]
            if role not in ("signal", "regions", "genes"):
                raise ValueError("Unknown file role")
            records.append(
                dict(
                    id=row["id"],
                    path=row["path"],
                    name=row["name"],
                    role=role,
                    source=row["source"],
                    group=str(item.get("group", row["grp"]) or "")[:80] or "group",
                    label=str(item.get("label") or row["label"])[:80],
                    color=str(item.get("color") or row["color"] or "")[:32],
                    cmap=str(item.get("cmap") or row["cmap"] or "")[:200],
                    track=row.get("track") or "",
                    info=row.get("info") or {},
                )
            )
        return records

    def upload_usage(self, project_id):
        row = self.db.one(
            "SELECT COALESCE(SUM(size), 0) AS n FROM files WHERE project_id = ? AND source = 'upload'",
            project_id,
        )
        return int(row["n"])

    def upload_dir(self, project_id):
        folder = self.uploads / project_id / new_id()
        folder.mkdir(parents=True)
        return folder

    def delete_file(self, row):
        self.db.run("DELETE FROM files WHERE id = ?", row["id"])
        folder = Path(row["path"]).parent
        if row["source"] == "upload" and folder.parent.parent == self.uploads:
            shutil.rmtree(folder, ignore_errors=True)


def _transcripts_gtf(bed: Path, gtf: Path) -> Path:
    """SeqPlots ships transcripts as BED; derive a minimal GTF for gene tracks."""
    if not bed.is_file() or (gtf.is_file() and gtf.stat().st_mtime >= bed.stat().st_mtime):
        return gtf
    from ..io.bed import read_bed

    rows = []
    for feature in read_bed(bed):
        interval = feature.interval
        name = feature.name or "tx"
        gene = name.rsplit(".", 1)[0]
        rows.append(
            f"{interval.chrom}\tseqplots\texon\t{interval.start + 1}\t{interval.end}\t.\t"
            f'{interval.strand or "."}\t.\tgene_id "{gene}"; transcript_id "{name}"; gene_name "{gene}";'
        )
    try:
        gtf.write_text("\n".join(rows) + "\n", encoding="utf-8")
    except OSError:  # read-only install (e.g. container): the shipped GTF is used as is
        pass
    return gtf
