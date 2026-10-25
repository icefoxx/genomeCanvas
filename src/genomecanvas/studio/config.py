"""Studio configuration: TOML file, then GC_STUDIO_* environment variables, then CLI flags."""

from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path


@dataclass
class Config:
    # Network
    host: str = "127.0.0.1"
    port: int = 8765
    public_url: str = (
        ""  # e.g. https://figures.core.example.edu; used for invite links and Origin checks
    )
    secure_cookies: bool = False  # set true behind HTTPS
    # Authentication: "none" (single local user), "local" (accounts in the database),
    # "proxy" (trusted reverse proxy such as oauth2-proxy passes the user's e-mail)
    auth: str = "none"
    proxy_user_header: str = "X-Auth-Request-Email"
    proxy_name_header: str = "X-Auth-Request-User"
    trusted_proxies: list[str] = field(default_factory=lambda: ["127.0.0.1", "::1"])
    proxy_auto_create: bool = True
    session_hours: int = 12
    # Storage
    data_root: str = (
        ""  # project deliverables, mounted read-only; project data_dir is relative to it
    )
    shared_root: str = ""  # shared annotations: <shared_root>/<genome>/*.gtf
    workspace: str = ""  # writable: database, uploads, results, cache
    static_dir: str = ""
    demo_project: bool = True  # bundled public example data available to every user
    # Jobs
    workers: int = 2
    job_timeout: int = 300
    worker_memory_mb: int = 4096
    max_jobs_per_user: int = 3
    max_queue: int = 64
    # Limits and retention
    max_upload_mb: int = 2048
    default_quota_gb: float = 20
    job_retention_days: int = 14
    upload_retention_days: int = 30
    cache_max_gb: float = 20

    @property
    def workspace_path(self) -> Path:
        return Path(self.workspace).expanduser().resolve()

    @property
    def db_path(self) -> Path:
        return self.workspace_path / "studio.db"


def _coerce(kind, value):
    if kind is bool or kind == "bool":
        return str(value).lower() in {"1", "true", "yes", "on"}
    if kind is int or kind == "int":
        return int(value)
    if kind is float or kind == "float":
        return float(value)
    if "list" in str(kind):
        return (
            [item.strip() for item in str(value).split(",") if item.strip()]
            if isinstance(value, str)
            else list(value)
        )
    return str(value)


def load_config(path: str | Path | None = None, overrides: dict | None = None) -> Config:
    config = Config()
    values: dict = {}
    if path:
        with open(path, "rb") as handle:
            data = tomllib.load(handle)
        values.update(data.get("studio", data))
    for item in fields(Config):
        env = os.environ.get(f"GC_STUDIO_{item.name.upper()}")
        if env is not None:
            values[item.name] = env
    values.update({k: v for k, v in (overrides or {}).items() if v is not None})
    known = {item.name: item.type for item in fields(Config)}
    for key, value in values.items():
        if key not in known:
            raise ValueError(f"Unknown studio setting '{key}'")
        setattr(config, key, _coerce(known[key], value))
    if config.auth not in {"none", "local", "proxy"}:
        raise ValueError("auth must be none, local or proxy")
    if not config.workspace:
        import tempfile

        config.workspace = str(Path(tempfile.gettempdir()) / "genomecanvas-studio")
    repo = Path(__file__).resolve().parents[3]
    if not config.static_dir:
        config.static_dir = str(repo / "web")
    if (
        config.auth != "none"
        and config.host not in {"127.0.0.1", "localhost", "::1"}
        and not config.secure_cookies
    ):
        print(
            "warning: serving accounts without secure cookies; put the studio behind HTTPS",
            flush=True,
        )
    return config


def example_config() -> str:
    return """# GenomeCanvas Studio configuration
[studio]
host = "127.0.0.1"          # keep on localhost behind nginx
port = 8765
public_url = "https://figures.core.example.edu"
secure_cookies = true
auth = "local"              # none | local | proxy
# auth = "proxy"            # with oauth2-proxy / Shibboleth in front:
# proxy_user_header = "X-Auth-Request-Email"
# trusted_proxies = ["127.0.0.1"]

data_root = "/data/core/deliverables"   # read-only; projects point at sub-folders
shared_root = "/data/core/annotations"  # <genome>/genes.gtf, cytoBand.txt ...
workspace = "/srv/genomecanvas"         # database, uploads, results, cache
demo_project = true

workers = 4
job_timeout = 300
worker_memory_mb = 4096
max_jobs_per_user = 3

max_upload_mb = 2048
default_quota_gb = 20
job_retention_days = 14
upload_retention_days = 30
cache_max_gb = 20
"""
