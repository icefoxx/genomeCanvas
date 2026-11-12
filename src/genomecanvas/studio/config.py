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
    # Chat assistant: the [agent] table (see genomecanvas.studio.agent.settings); empty = off
    agent: dict = field(default_factory=dict)

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
        values["agent"] = data.get("agent", values.get("agent")) or {}
    for item in fields(Config):
        env = os.environ.get(f"GC_STUDIO_{item.name.upper()}")
        if env is not None:
            values[item.name] = env
    values.update({k: v for k, v in (overrides or {}).items() if v is not None})
    known = {item.name: item.type for item in fields(Config)}
    agent = values.pop("agent", None)
    for key, value in values.items():
        if key not in known:
            raise ValueError(f"Unknown studio setting '{key}'")
        setattr(config, key, _coerce(known[key], value))
    if isinstance(agent, dict):
        from .agent.settings import agent_settings

        agent_settings(agent)  # fail at start-up, not on the first chat message
        config.agent = agent
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

# Chat assistant (pip install -e ".[agent]" for Claude). Leave out to turn it off.
# API keys come from environment variables, never from this file.
# [agent]
# max_steps = 12                       # tool calls per message
# daily_tokens_per_user = 2000000      # 0 = no limit
# local_only_projects = ["P2026-017"]  # these may only use providers with local = true
#
# [[agent.providers]]
# name = "claude"
# kind = "anthropic"
# models = ["claude-opus-5-5", "claude-sonnet-5-5"]
# api_key_env = "ANTHROPIC_API_KEY"
#
# [[agent.providers]]
# name = "openai"
# kind = "openai"
# models = ["<model id from your OpenAI account>"]
# api_key_env = "OPENAI_API_KEY"
#
# [[agent.providers]]
# name = "lab-gpu"                     # vLLM / Ollama / LM Studio on site
# kind = "openai"
# base_url = "http://gpu01:8000/v1"
# models = ["<served model name>"]
# local = true
# vision = false                       # text-only models get the methods text, not the image
#
# [[agent.providers]]
# name = "opencode"                    # OpenCode Zen gateway, key in OPENCODE_API_KEY
# kind = "opencode"                    # plan = "go" for an OpenCode Go subscription
# models = ["kimi-k3", "glm-5.3"]      # paid Zen text models; *-free ones only work in OpenCode
# ca_bundle = "system"                 # behind HTTPS inspection (e.g. Cloudflare Gateway)
#
# [[agent.providers]]
# name = "opencode-claude"             # Claude models through Zen's Anthropic endpoint
# kind = "anthropic"
# base_url = "https://opencode.ai/inference/anthropic"
# api_key_env = "OPENCODE_API_KEY"
# auth = "bearer"                      # send the key as Authorization: Bearer
# models = ["claude-opus-5-5"]
"""
