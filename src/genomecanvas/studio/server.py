"""GenomeCanvas Studio HTTP server (stdlib; run behind nginx/HTTPS for multi-user use)."""

from __future__ import annotations

import ipaddress
import json
import mimetypes
import re
import threading
import time
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

from .config import Config
from .db import Database, check_password, new_id
from .engine import CATALOG, SUFFIXES, gene_search, inspect_file, meta
from .jobs import JobRunner, friendly
from .projects import SHARED, Projects

STATIC = {
    "/": "index.html",
    "/index.html": "index.html",
    "/app.js": "app.js",
    "/styles.css": "styles.css",
    "/login": "login.html",
    "/invite": "login.html",
    "/login.js": "login.js",
    "/admin": "admin.html",
    "/admin.js": "admin.js",
}
COOKIE = "gc_session"
SECURITY_HEADERS = {
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "DENY",
    "Referrer-Policy": "same-origin",
    "Content-Security-Policy": (
        "default-src 'self'; img-src 'self' data: blob:; style-src 'self' 'unsafe-inline'; "
        "script-src 'self'; connect-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'"
    ),
}


class HttpError(Exception):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


class Studio:
    """Application state shared by request handler threads."""

    def __init__(self, config: Config):
        self.config = config
        config.workspace_path.mkdir(parents=True, exist_ok=True)
        self.db = Database(config.db_path)
        self.projects = Projects(self.db, config)
        self.jobs = JobRunner(self.db, config)
        self.static = Path(config.static_dir)
        self.login_attempts: dict[str, list[float]] = {}
        self.attempts_lock = threading.Lock()
        self.meta = None

    def start(self):
        self.projects.ensure_demos()
        self.projects.scan_shared()
        if self.config.auth == "none" and not self.db.user_by_email("local"):
            self.db.create_user("local", name="Local user", is_admin=True)
        self.jobs.start()
        threading.Thread(target=self._scan_all, daemon=True).start()

    def _scan_all(self):
        for project in self.db.all(
            "SELECT id FROM projects WHERE status = 'active' AND data_dir != ''"
        ):
            try:
                self.projects.scan(project["id"])
            except Exception as exc:  # noqa: BLE001
                print(f"scan {project['id']}: {exc}", flush=True)

    LIMITS = {"account": 10, "ip": 60}  # failed attempts per 10 minutes

    def _recent(self, key, now):
        recent = [t for t in self.login_attempts.get(key, []) if now - t < 600]
        self.login_attempts[key] = recent
        return recent

    def check_throttle(self, ip, account):
        """Refuse when too many recent *failures* for this account+IP or for the IP overall."""
        now = time.time()
        with self.attempts_lock:
            too_many = (
                len(self._recent(f"{ip}|{account}", now)) >= self.LIMITS["account"]
                or len(self._recent(ip, now)) >= self.LIMITS["ip"]
            )
        if too_many:
            raise HttpError(429, "Too many failed sign-in attempts; wait a few minutes")

    def record_failure(self, ip, account):
        now = time.time()
        with self.attempts_lock:
            for key in (f"{ip}|{account}", ip):
                self._recent(key, now).append(now)

    def trusted(self, address):
        """Is a connection coming from a configured reverse proxy (IP or CIDR)?"""
        try:
            addr = ipaddress.ip_address(address)
        except ValueError:
            return False
        for entry in self.config.trusted_proxies:
            try:
                if addr in ipaddress.ip_network(entry, strict=False):
                    return True
            except ValueError:
                continue
        return False


def public_file(row, role):
    shared = row["project_id"] == SHARED
    return dict(
        id=row["id"],
        name=row["name"],
        role=row["role"],
        label=row["label"],
        group=row["grp"] or "",
        color=row["color"] or "",
        cmap=row["cmap"] or "",
        source=row["source"],
        info=row["info"] or {},
        shared=shared,
        deletable=row["source"] == "upload" and role in ("admin", "editor"),
    )


def public_user(user):
    return {
        k: user[k]
        for k in ("id", "email", "name", "is_admin", "active", "created", "last_login")
        if k in user
    }


def make_handler(app: Studio):
    config = app.config

    class Handler(BaseHTTPRequestHandler):
        server_version = "GenomeCanvasStudio"
        sys_version = ""

        def log_message(self, format, *args):  # noqa: A002
            if "/api/jobs/" not in (args[0] if args else ""):
                super().log_message(format, *args)

        # ---------------------------------------------------- plumbing
        @property
        def ip(self):
            """Real client address. Behind a trusted proxy use the right-most X-Forwarded-For
            entry: it is the one our proxy appended, whereas earlier entries are client-supplied."""
            peer = self.client_address[0]
            forwarded = self.headers.get("X-Forwarded-For")
            if forwarded and app.trusted(peer):
                return forwarded.split(",")[-1].strip() or peer
            return peer

        def send(self, data, status=200, content_type="application/json", headers=None):
            body = data if isinstance(data, bytes) else json.dumps(data, default=str).encode()
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            for key, value in {**SECURITY_HEADERS, **(headers or {})}.items():
                self.send_header(key, value)
            self.end_headers()
            self.wfile.write(body)

        def redirect(self, location):
            self.send_response(303)
            self.send_header("Location", location)
            self.send_header("Content-Length", "0")
            self.end_headers()

        def json_body(self, limit=4 * 1024 * 1024):
            length = int(self.headers.get("Content-Length") or 0)
            if not 0 < length <= limit:
                raise HttpError(400, "Request body missing or too large")
            try:
                body = json.loads(self.rfile.read(length))
            except json.JSONDecodeError as exc:
                raise HttpError(400, f"Invalid JSON: {exc}") from exc
            if not isinstance(body, dict):
                raise HttpError(400, "Expected a JSON object")
            return body

        def cookie_token(self):
            raw = self.headers.get("Cookie")
            if not raw:
                return None
            morsel = SimpleCookie(raw).get(COOKIE)
            return morsel.value if morsel else None

        def session_cookie(self, token, max_age):
            parts = [
                f"{COOKIE}={token}",
                "Path=/",
                "HttpOnly",
                "SameSite=Lax",
                f"Max-Age={max_age}",
            ]
            if config.secure_cookies:
                parts.append("Secure")
            return {"Set-Cookie": "; ".join(parts)}

        def user(self, required=True):
            user = None
            if config.auth == "none":
                user = app.db.user_by_email("local")
            elif config.auth == "proxy":
                if app.trusted(self.client_address[0]):
                    email = (self.headers.get(config.proxy_user_header) or "").strip().lower()
                    if email:
                        user = app.db.user_by_email(email)
                        if user is None and config.proxy_auto_create:
                            app.db.create_user(
                                email, name=self.headers.get(config.proxy_name_header) or None
                            )
                            user = app.db.user_by_email(email)
                            app.db.audit(user, "user.auto_created", ip=self.ip)
                        if user and not user["active"]:
                            user = None
            else:
                user = app.db.session_user(self.cookie_token())
            if required and user is None:
                raise HttpError(401, "Please sign in")
            return user

        def guard_write(self):
            """CSRF defence: custom header (forces CORS preflight) plus same-origin check."""
            if self.headers.get("X-GenomeCanvas") != "studio":
                raise HttpError(403, "Missing studio request header")
            origin = self.headers.get("Origin")
            if origin:
                allowed = {
                    f"http://{self.headers.get('Host')}",
                    f"https://{self.headers.get('Host')}",
                }
                if config.public_url:
                    allowed.add(config.public_url.rstrip("/"))
                if origin not in allowed:
                    raise HttpError(403, "Cross-origin requests are not allowed")

        def project_access(self, user, project_id, need="viewer"):
            role, project = app.projects.role(user, project_id)
            if project is None or role is None:
                raise HttpError(404, "Project not found or access has ended")
            if need == "editor" and role not in ("editor", "admin"):
                raise HttpError(403, "Your access to this project is view-only")
            if need == "admin" and role != "admin":
                raise HttpError(403, "Administrators only")
            return role, project

        def require_admin(self, user):
            if not user["is_admin"]:
                raise HttpError(403, "Administrators only")

        def dispatch(self, method):
            parsed = urlparse(self.path)
            path = unquote(parsed.path)
            query = {k: v[0] for k, v in parse_qs(parsed.query).items()}
            try:
                if method != "GET":
                    self.guard_write()
                return self.route(method, path, query)
            except HttpError as exc:
                if (
                    exc.status == 401
                    and method == "GET"
                    and not path.startswith(("/api/", "/results/"))
                ):
                    return self.redirect("/login")
                return self.send({"error": str(exc)}, exc.status)
            except (ValueError, KeyError, TypeError, FileNotFoundError) as exc:
                return self.send({"error": friendly(exc)}, 400)

        def do_GET(self):  # noqa: N802
            self.dispatch("GET")

        def do_POST(self):  # noqa: N802
            self.dispatch("POST")

        def do_PATCH(self):  # noqa: N802
            self.dispatch("PATCH")

        def do_DELETE(self):  # noqa: N802
            self.dispatch("DELETE")

        # ---------------------------------------------------- routing
        def route(self, method, path, query):  # noqa: C901 - flat routing table reads best here
            if method == "GET" and path in STATIC:
                name = STATIC[path]
                if name in ("index.html", "admin.html") and config.auth != "none":
                    user = self.user()
                    if name == "admin.html":
                        self.require_admin(user)
                if name == "login.html" and config.auth != "local":
                    return self.redirect("/")
                file = app.static / name
                if not file.is_file():
                    raise HttpError(404, "Not found")
                return self.send(
                    file.read_bytes(), content_type=mimetypes.guess_type(name)[0] or "text/plain"
                )
            if path == "/api/login" and method == "POST":
                return self.login()
            if path == "/api/invite" and method == "POST":
                return self.accept_invite()
            if path == "/api/logout" and method == "POST":
                token = self.cookie_token()
                if token:
                    from .db import token_hash

                    app.db.run("DELETE FROM sessions WHERE token_hash = ?", token_hash(token))
                return self.send({"ok": True}, headers=self.session_cookie("", 0))
            if path == "/api/health":
                return self.send({"status": "ok"})
            if path == "/vendor/plotly.min.js" and method == "GET":
                return self.vendor_plotly()

            user = self.user()
            parts = [p for p in path.split("/") if p]
            if path == "/api/me":
                return self.send(
                    dict(
                        user=public_user(user),
                        auth=config.auth,
                        projects=[self.project_summary(p) for p in app.projects.accessible(user)],
                        limits=dict(
                            max_upload=config.max_upload_mb * 1024 * 1024,
                            job_timeout=config.job_timeout,
                        ),
                    )
                )
            if path == "/api/meta":
                if app.meta is None:
                    app.meta = meta()
                return self.send(dict(app.meta, max_upload=config.max_upload_mb * 1024 * 1024))
            if parts[:2] == ["api", "projects"] and len(parts) >= 3:
                return self.project_routes(method, user, parts[2], parts[3:], query)
            if parts[:2] == ["api", "jobs"] and len(parts) == 3 and method == "GET":
                return self.job_status(user, parts[2])
            if parts[:2] == ["api", "figures"] and len(parts) == 3:
                return self.figure_item(method, user, parts[2])
            if parts[:2] == ["api", "templates"]:
                return self.template_routes(method, user, parts[2:])
            if parts[:2] == ["api", "admin"]:
                self.require_admin(user)
                return self.admin_routes(method, user, parts[2:], query)
            if parts[:1] == ["results"] and len(parts) == 3 and method == "GET":
                return self.result_file(user, parts[1], parts[2], query)
            raise HttpError(404, "Not found")

        # ---------------------------------------------------- auth endpoints
        def login(self):
            if config.auth != "local":
                raise HttpError(400, "Password sign-in is not enabled")
            body = self.json_body(64 * 1024)
            account = str(body.get("email", "")).strip().lower()[:200]
            app.check_throttle(self.ip, account)
            user = app.db.user_by_email(account)
            if (
                not user
                or not user["active"]
                or not check_password(str(body.get("password", "")), user["password_hash"])
            ):
                app.record_failure(self.ip, account)
                app.db.audit(user, "login.failed", detail=account[:120], ip=self.ip)
                raise HttpError(401, "E-mail or password is incorrect")
            token = app.db.create_session(user["id"], config.session_hours)
            app.db.audit(user, "login", ip=self.ip)
            return self.send(
                {"ok": True}, headers=self.session_cookie(token, config.session_hours * 3600)
            )

        def accept_invite(self):
            if config.auth != "local":
                raise HttpError(400, "Invitations are only used with password sign-in")
            body = self.json_body(64 * 1024)
            app.check_throttle(self.ip, "invite")
            try:
                user = app.db.redeem_invite(
                    str(body.get("token", "")), str(body.get("password", ""))
                )
            except ValueError:
                app.record_failure(self.ip, "invite")
                raise
            token = app.db.create_session(user["id"], config.session_hours)
            app.db.audit(user, "invite.accepted", ip=self.ip)
            return self.send(
                {"ok": True}, headers=self.session_cookie(token, config.session_hours * 3600)
            )

        # ---------------------------------------------------- projects
        def project_summary(self, project):
            return dict(
                id=project["id"],
                name=project["name"],
                description=project.get("description") or "",
                genome=project.get("genome") or "",
                role=project.get("role"),
                demo=project.get("demo", False),
                expires=project.get("expires"),
                status=project.get("status"),
            )

        def project_routes(self, method, user, project_id, rest, query):
            need = "viewer" if method == "GET" else "editor"
            role, project = self.project_access(user, project_id, need)
            files = lambda: app.projects.files(project)  # noqa: E731
            if not rest and method == "GET":
                templates = app.db.all(
                    "SELECT * FROM templates WHERE project_id IS NULL OR project_id = ? ORDER BY project_id IS NULL, name",
                    project_id,
                )
                return self.send(
                    dict(
                        project=dict(
                            self.project_summary(
                                {**project, "role": role, "demo": app.projects.is_demo(project_id)}
                            ),
                            defaults=project.get("defaults") or {},
                            quota=project.get("quota_bytes"),
                            used=app.projects.upload_usage(project_id),
                        ),
                        role=role,
                        files=[public_file(f, role) for f in files()],
                        templates=[self.public_template(t) for t in templates],
                    )
                )
            if rest == ["files"] and method == "GET":
                return self.send({"files": [public_file(f, role) for f in files()]})
            if rest == ["uploads"] and method == "POST":
                return self.upload(user, project)
            if len(rest) == 2 and rest[0] == "files" and method == "DELETE":
                row = app.db.one(
                    "SELECT * FROM files WHERE id = ? AND project_id = ?", rest[1], project_id
                )
                if not row or row["source"] != "upload":
                    raise HttpError(404, "Only uploaded files can be deleted")
                app.projects.delete_file(row)
                app.db.audit(user, "file.deleted", project_id, row["name"], self.ip)
                return self.send({"deleted": row["id"]})
            if rest == ["genes"] and method == "GET":
                wanted = set((query.get("files") or "").split(","))
                paths = [f["path"] for f in files() if f["role"] == "genes" and f["id"] in wanted]
                return self.send({"hits": gene_search(paths, query.get("q", ""))})
            if rest == ["jobs"] and method == "POST":
                return self.submit_job(user, project)
            if rest == ["figures"] and method == "GET":
                rows = app.db.all(
                    "SELECT f.*, u.email AS author, j.status AS job_status FROM figures f "
                    "LEFT JOIN users u ON u.id = f.created_by LEFT JOIN jobs j ON j.id = f.job_id "
                    "WHERE f.project_id = ? ORDER BY f.updated DESC",
                    project_id,
                )
                return self.send({"figures": [self.public_figure(r) for r in rows]})
            if rest == ["figures"] and method == "POST":
                body = self.json_body()
                name = str(body.get("name") or "").strip()[:120]
                if not name:
                    raise HttpError(400, "Give the figure a name")
                job_id = body.get("job_id")
                if job_id and not app.db.one(
                    "SELECT id FROM jobs WHERE id = ? AND project_id = ? AND status = 'complete'",
                    job_id,
                    project_id,
                ):
                    raise HttpError(
                        400, "That figure's result is no longer available; plot again first"
                    )
                now = time.time()
                figure = dict(
                    id=new_id(),
                    project_id=project_id,
                    name=name,
                    notes=str(body.get("notes") or "")[:4000],
                    config=body.get("config") or {},
                    job_id=job_id,
                    created_by=user["id"],
                    created=now,
                    updated=now,
                )
                app.db.insert("figures", **figure)
                app.db.audit(user, "figure.saved", project_id, name, self.ip)
                return self.send(
                    {
                        "figure": self.public_figure(
                            {**figure, "author": user["email"], "job_status": "complete"}
                        )
                    },
                    201,
                )
            raise HttpError(404, "Not found")

        def upload(self, user, project):
            length = int(self.headers.get("Content-Length") or 0)
            if not 0 < length <= config.max_upload_mb * 1024 * 1024:
                raise HttpError(413, f"Uploads are limited to {config.max_upload_mb} MB per file")
            quota = project.get("quota_bytes") or 0
            used = app.projects.upload_usage(project["id"])
            if quota and used + length > quota:
                raise HttpError(
                    413,
                    f"Project upload quota ({quota / 1024**3:.1f} GB) would be exceeded; delete old uploads first",
                )
            name = Path(unquote(self.headers.get("X-Filename", ""))).name
            suffix = Path(name).suffix.lower()
            if suffix == ".gz":
                raise HttpError(
                    400,
                    "Please decompress .gz files first (BED, bedGraph and GTF are read as plain text)",
                )
            if suffix not in SUFFIXES or not re.fullmatch(r"[\w.+ ()-]{1,200}", name):
                raise HttpError(
                    400,
                    "Supported files: BigWig (.bw), bedGraph, BED / narrowPeak / broadPeak, GTF",
                )
            folder = app.projects.upload_dir(project["id"])
            target = folder / name
            remaining = length
            with target.open("wb") as handle:
                while remaining:
                    chunk = self.rfile.read(min(1 << 20, remaining))
                    if not chunk:
                        break
                    handle.write(chunk)
                    remaining -= len(chunk)
            role = SUFFIXES[suffix]
            info = (
                inspect_file(target, role) if not remaining else {"error": "Upload was interrupted"}
            )
            if info.get("error"):
                import shutil

                shutil.rmtree(folder, ignore_errors=True)
                raise HttpError(400, f"{name}: {info['error']}")
            row = dict(
                id=new_id(),
                project_id=project["id"],
                path=str(target),
                name=name,
                role=role,
                label=name.split(".")[0],
                grp="upload",
                color="",
                cmap="",
                track="",
                source="upload",
                size=length,
                mtime=target.stat().st_mtime,
                sort=0,
                info=info,
                uploaded_by=user["id"],
                created=time.time(),
            )
            app.db.insert("files", **row)
            app.db.audit(
                user, "file.uploaded", project["id"], dict(name=name, size=length), self.ip
            )
            return self.send(public_file({**row, "info": info}, "editor"), 201)

        # ---------------------------------------------------- jobs & results
        def submit_job(self, user, project):
            body = self.json_body()
            records = app.projects.resolve(project, body.get("files"))
            settings = body.get("settings") if isinstance(body.get("settings"), dict) else {}
            if settings.get("rel_blacklist"):
                # The blacklist is authorised like any other input, then tagged so it is not plotted.
                (black,) = app.projects.resolve(
                    project, [{"id": settings["rel_blacklist"], "role": "regions"}]
                )
                records.append({**black, "role": "blacklist"})
            active = app.db.one(
                "SELECT COUNT(*) AS n FROM jobs WHERE user_id = ? AND status IN ('queued','running')",
                user["id"],
            )["n"]
            if active >= config.max_jobs_per_user:
                raise HttpError(
                    429, "You already have figures being drawn; wait for them to finish"
                )
            if (
                app.db.one("SELECT COUNT(*) AS n FROM jobs WHERE status IN ('queued','running')")[
                    "n"
                ]
                >= config.max_queue
            ):
                raise HttpError(429, "The studio is busy; try again in a minute")
            payload = {
                k: body.get(k)
                for k in ("category", "plot", "settings", "title", "palette", "advanced")
            }
            payload["files"] = [
                {k: r[k] for k in ("id", "role", "group", "label", "color", "cmap")}
                for r in records
            ]
            job_id = new_id()
            app.db.insert(
                "jobs",
                id=job_id,
                project_id=project["id"],
                user_id=user["id"],
                status="queued",
                payload=payload,
                created=time.time(),
            )
            app.jobs.submit(job_id, project["id"], payload, records)
            category = payload.get("category")
            plot = payload.get("plot")
            try:
                name = CATALOG[int(category)][1][int(plot)]
            except (TypeError, ValueError, IndexError):
                name = "figure"
            app.db.audit(user, "job.submitted", project["id"], name, self.ip)
            return self.send({"id": job_id, "status": "queued"}, 202)

        def job_status(self, user, job_id):
            job = app.db.one("SELECT * FROM jobs WHERE id = ?", job_id)
            if not job:
                raise HttpError(404, "Job not found")
            self.project_access(user, job["project_id"])
            out = dict(id=job_id, status=job["status"], error=job["error"], elapsed=job["elapsed"])
            if job["status"] == "queued":
                out["position"] = app.jobs.position(job_id)
            if job["status"] == "complete":
                directory = app.jobs.job_dir(job["project_id"], job_id)
                result = job["result"] or {}
                out.update(
                    warnings=result.get("warnings", []),
                    counts=result.get("counts"),
                    files=[
                        {
                            "name": p.name,
                            "url": f"/results/{job_id}/{p.name}",
                            "size": p.stat().st_size,
                        }
                        for p in sorted(directory.iterdir())
                        if p.is_file()
                    ]
                    if directory.is_dir()
                    else [],
                    preview=f"/results/{job_id}/figure.svg",
                    preview_png=f"/results/{job_id}/figure.png",
                )
            return self.send(out)

        def result_file(self, user, job_id, name, query):
            if not re.fullmatch(r"[a-f0-9]{32}", job_id) or Path(name).name != name:
                raise HttpError(404, "Invalid artifact")
            job = app.db.one("SELECT project_id FROM jobs WHERE id = ?", job_id)
            if not job:
                raise HttpError(404, "This result has expired")
            self.project_access(user, job["project_id"])
            file = app.jobs.job_dir(job["project_id"], job_id) / name
            if not file.is_file():
                raise HttpError(404, "Artifact not found")
            content_type = mimetypes.guess_type(name)[0] or "application/octet-stream"
            if file.suffix in (".tsv", ".txt", ".bed", ".yaml"):
                content_type = "text/plain; charset=utf-8"
            headers = {}
            if file.suffix in (".html", ".htm"):
                # Result HTML runs its own scripts: never render it inside the studio origin.
                query = {**query, "download": query.get("download") or name}
            if query.get("download"):
                safe = re.sub(r"[^\w.-]+", "_", query["download"])[:120] or name
                headers["Content-Disposition"] = f'attachment; filename="{safe}"'
            if file.suffix == ".svg":
                headers["Content-Security-Policy"] = "default-src 'none'; style-src 'unsafe-inline'"
            return self.send(file.read_bytes(), content_type=content_type, headers=headers)

        # ---------------------------------------------------- figures & templates
        def public_figure(self, row):
            available = row.get("job_status") == "complete"
            return dict(
                id=row["id"],
                name=row["name"],
                notes=row.get("notes") or "",
                config=row.get("config") or {},
                job_id=row.get("job_id") if available else None,
                author=row.get("author"),
                created=row.get("created"),
                updated=row.get("updated"),
                thumb=f"/results/{row['job_id']}/figure.png"
                if available and row.get("job_id")
                else None,
            )

        def figure_item(self, method, user, figure_id):
            figure = app.db.one("SELECT * FROM figures WHERE id = ?", figure_id)
            if not figure:
                raise HttpError(404, "Figure not found")
            role, _ = self.project_access(
                user, figure["project_id"], "viewer" if method == "GET" else "editor"
            )
            if method == "PATCH":
                body = self.json_body()
                updates = {k: str(body[k])[:4000] for k in ("name", "notes") if k in body}
                if "config" in body and isinstance(body["config"], dict):
                    updates["config"] = body["config"]
                if body.get("job_id"):
                    if not app.db.one(
                        "SELECT id FROM jobs WHERE id = ? AND project_id = ? AND status = 'complete'",
                        body["job_id"],
                        figure["project_id"],
                    ):
                        raise HttpError(400, "Plot the figure before updating it")
                    updates["job_id"] = body["job_id"]
                app.db.update("figures", "id", figure_id, updated=time.time(), **updates)
                app.db.audit(user, "figure.updated", figure["project_id"], figure["name"], self.ip)
                return self.send({"ok": True})
            if method == "DELETE":
                if role != "admin" and figure["created_by"] != user["id"]:
                    raise HttpError(
                        403, "Only the author or an administrator can delete this figure"
                    )
                app.db.run("DELETE FROM figures WHERE id = ?", figure_id)
                app.db.audit(user, "figure.deleted", figure["project_id"], figure["name"], self.ip)
                return self.send({"deleted": figure_id})
            return self.send({"figure": figure})

        def public_template(self, row):
            return dict(
                id=row["id"],
                name=row["name"],
                description=row.get("description") or "",
                config=row.get("config") or {},
                scope="project" if row.get("project_id") else "global",
                managed=row.get("created_by") == "manifest",
            )

        def template_routes(self, method, user, rest):
            self.require_admin(user)
            if not rest and method == "POST":
                body = self.json_body()
                name = str(body.get("name") or "").strip()[:120]
                config_ = body.get("config")
                if not name or not isinstance(config_, dict):
                    raise HttpError(400, "Templates need a name and a configuration")
                project_id = body.get("project_id") or None
                if project_id and not app.db.one(
                    "SELECT id FROM projects WHERE id = ?", project_id
                ):
                    raise HttpError(404, "Unknown project")
                row = dict(
                    id=new_id(),
                    project_id=project_id,
                    name=name,
                    description=str(body.get("description") or "")[:500],
                    config={
                        k: config_.get(k) for k in ("category", "plot", "settings", "advanced")
                    },
                    created_by=user["id"],
                    created=time.time(),
                )
                app.db.insert("templates", **row)
                app.db.audit(user, "template.created", project_id, name, self.ip)
                return self.send({"template": self.public_template(row)}, 201)
            if len(rest) == 1 and method == "DELETE":
                app.db.run(
                    "DELETE FROM templates WHERE id = ? AND created_by != 'manifest'", rest[0]
                )
                app.db.audit(user, "template.deleted", None, rest[0], self.ip)
                return self.send({"deleted": rest[0]})
            raise HttpError(404, "Not found")

        # ---------------------------------------------------- admin
        def admin_routes(self, method, user, rest, query):  # noqa: C901
            db = app.db
            if rest == ["overview"] and method == "GET":
                users = [
                    public_user(u) for u in db.all("SELECT * FROM users ORDER BY created DESC")
                ]
                since = time.time() - 30 * 86400
                projects = []
                for p in db.all("SELECT * FROM projects ORDER BY created DESC"):
                    members = db.all(
                        "SELECT u.id, u.email, u.name, m.role FROM members m JOIN users u ON u.id = m.user_id "
                        "WHERE m.project_id = ? ORDER BY u.email",
                        p["id"],
                    )
                    usage = db.one(
                        "SELECT COUNT(*) AS jobs, COALESCE(SUM(bytes),0) AS result_bytes, "
                        "SUM(status = 'error') AS errors, MAX(created) AS last_job FROM jobs WHERE project_id = ? AND created > ?",
                        p["id"],
                        since,
                    )
                    files = db.one(
                        "SELECT SUM(source != 'upload') AS delivered, SUM(source = 'upload') AS uploads FROM files WHERE project_id = ?",
                        p["id"],
                    )
                    projects.append(
                        dict(
                            p,
                            members=members,
                            usage=usage,
                            files=files,
                            upload_bytes=app.projects.upload_usage(p["id"]),
                            figures=db.one(
                                "SELECT COUNT(*) AS n FROM figures WHERE project_id = ?", p["id"]
                            )["n"],
                            demo=p["id"] in ("demo", "demo-ce10"),
                        )
                    )
                return self.send(
                    dict(
                        users=users,
                        projects=projects,
                        auth=config.auth,
                        data_root=str(app.projects.data_root or ""),
                        shared_root=str(app.projects.shared_root or ""),
                        templates=[
                            self.public_template(t)
                            for t in db.all("SELECT * FROM templates ORDER BY created DESC")
                        ],
                    )
                )
            if rest == ["audit"] and method == "GET":
                limit = min(2000, int(query.get("limit", 300)))
                rows = db.all(
                    "SELECT a.*, u.email FROM audit a LEFT JOIN users u ON u.id = a.user_id ORDER BY a.id DESC LIMIT ?",
                    limit,
                )
                return self.send({"audit": rows})
            if rest == ["users"] and method == "POST":
                body = self.json_body()
                created = db.create_user(
                    str(body.get("email", "")),
                    name=body.get("name") or None,
                    is_admin=bool(body.get("is_admin")),
                )
                db.audit(user, "user.created", detail=created["email"], ip=self.ip)
                return self.send(
                    {"user": public_user(created), **self.invite_payload(created)}, 201
                )
            if len(rest) == 2 and rest[0] == "users" and method == "PATCH":
                body = self.json_body()
                updates = {}
                if "active" in body:
                    updates["active"] = int(bool(body["active"]))
                if "is_admin" in body:
                    if rest[1] == user["id"] and not body["is_admin"]:
                        raise HttpError(400, "You cannot remove your own administrator role")
                    updates["is_admin"] = int(bool(body["is_admin"]))
                if "name" in body:
                    updates["name"] = str(body["name"])[:120]
                db.update("users", "id", rest[1], **updates)
                if updates.get("active") == 0:
                    db.run("DELETE FROM sessions WHERE user_id = ?", rest[1])
                db.audit(user, "user.updated", detail=dict(id=rest[1], **updates), ip=self.ip)
                return self.send({"ok": True})
            if len(rest) == 3 and rest[0] == "users" and rest[2] == "invite" and method == "POST":
                target = db.one("SELECT * FROM users WHERE id = ?", rest[1])
                if not target:
                    raise HttpError(404, "User not found")
                db.audit(user, "user.invited", detail=target["email"], ip=self.ip)
                return self.send(self.invite_payload(target))
            if rest == ["projects"] and method == "POST":
                body = self.json_body()
                project = app.projects.create(
                    str(body.get("id", "")).strip(),
                    name=body.get("name"),
                    data_dir=str(body.get("data_dir") or "").strip(),
                    genome=str(body.get("genome") or "").strip(),
                    description=str(body.get("description") or ""),
                    expires=body.get("expires") or None,
                    quota_gb=float(body["quota_gb"])
                    if body.get("quota_gb") not in (None, "")
                    else None,
                )
                report = app.projects.scan(project["id"]) if project.get("data_dir") else {}
                db.audit(user, "project.created", project["id"], ip=self.ip)
                return self.send({"project": project, "scan": report}, 201)
            if len(rest) == 2 and rest[0] == "projects" and method == "PATCH":
                body = self.json_body()
                project = db.one("SELECT * FROM projects WHERE id = ?", rest[1])
                if not project:
                    raise HttpError(404, "Project not found")
                updates = {
                    k: body[k]
                    for k in ("name", "description", "genome", "data_dir", "expires")
                    if k in body
                }
                if "status" in body:
                    if body["status"] not in ("active", "archived"):
                        raise HttpError(400, "Status must be active or archived")
                    updates["status"] = body["status"]
                if body.get("quota_gb") not in (None, ""):
                    updates["quota_bytes"] = int(float(body["quota_gb"]) * 1024**3)
                if updates.get("data_dir"):
                    app.projects.project_dir(dict(data_dir=updates["data_dir"]))
                if "expires" in updates and not updates["expires"]:
                    updates["expires"] = None
                db.update("projects", "id", rest[1], **updates)
                db.audit(user, "project.updated", rest[1], updates, self.ip)
                return self.send({"ok": True})
            if len(rest) == 3 and rest[0] == "projects" and rest[2] == "scan" and method == "POST":
                report = app.projects.scan(rest[1])
                db.audit(user, "project.scanned", rest[1], report, self.ip)
                return self.send({"scan": report})
            if (
                len(rest) == 3
                and rest[0] == "projects"
                and rest[2] == "members"
                and method == "POST"
            ):
                body = self.json_body()
                if not db.one("SELECT id FROM projects WHERE id = ?", rest[1]):
                    raise HttpError(404, "Project not found")
                email = str(body.get("email", "")).strip().lower()
                target = db.user_by_email(email)
                invite = {}
                if target is None:
                    if not body.get("role"):
                        raise HttpError(404, "User not found")
                    target = db.create_user(email)
                    invite = self.invite_payload(target)
                    db.audit(user, "user.created", detail=email, ip=self.ip)
                app.projects.set_member(rest[1], target["id"], body.get("role") or None)
                db.audit(
                    user, "member.set", rest[1], dict(email=email, role=body.get("role")), self.ip
                )
                return self.send({"ok": True, **invite})
            if rest == ["cleanup"] and method == "POST":
                report = app.jobs.cleanup()
                db.audit(user, "cleanup", detail=report, ip=self.ip)
                return self.send({"cleanup": report})
            if rest == ["rescan-shared"] and method == "POST":
                return self.send({"scan": app.projects.scan_shared()})
            raise HttpError(404, "Not found")

        def vendor_plotly(self):
            """Plotly's own bundle from the installed Python package (no CDN under our CSP)."""
            try:
                import plotly

                bundle = Path(plotly.__file__).parent / "package_data" / "plotly.min.js"
                body = bundle.read_bytes()
            except (ImportError, OSError) as exc:
                raise HttpError(
                    404, "Interactive plots need the 'plotly' Python package on the server"
                ) from exc
            self.send_response(200)
            self.send_header("Content-Type", "text/javascript; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "public, max-age=604800")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self.wfile.write(body)

        def invite_payload(self, target):
            if config.auth != "local":
                return {
                    "invite_url": None,
                    "note": f"Users sign in through {config.auth} authentication",
                }
            token = app.db.create_invite(target["id"])
            base = config.public_url.rstrip("/") or f"http://{self.headers.get('Host')}"
            return {"invite_url": f"{base}/invite?token={token}", "expires_days": 7}

    return Handler


def serve(config: Config):
    app = Studio(config)
    app.start()
    httpd = ThreadingHTTPServer((config.host, config.port), make_handler(app))
    httpd.daemon_threads = True
    print(
        f"GenomeCanvas Studio ({config.auth} auth) at http://{config.host}:{config.port}",
        flush=True,
    )
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        app.jobs.stop()
        httpd.server_close()
