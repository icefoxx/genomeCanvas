"""Chat-agent settings: the `[agent]` table of studio.toml.

    [agent]
    max_steps = 12                     # tool calls per user message
    daily_tokens_per_user = 2000000    # 0 = no limit
    local_only_projects = ["P2026-017"] # these projects may only use providers with local = true

    [[agent.providers]]
    name = "claude"
    kind = "anthropic"                 # anthropic | openai (any OpenAI-compatible endpoint)
    models = ["claude-opus-5-5", "claude-sonnet-5-5"]
    api_key_env = "ANTHROPIC_API_KEY"

    [[agent.providers]]
    name = "lab-gpu"
    kind = "openai"
    base_url = "http://gpu01:8000/v1"  # vLLM / Ollama / LM Studio
    models = ["qwen3-72b"]
    local = true
    vision = false

    [[agent.providers]]
    name = "opencode"
    kind = "opencode"                  # OpenCode Zen (pay as you go); plan = "go" for OpenCode Go
    models = ["kimi-k3", "glm-5.3"]    # any Zen text model; *-free models only work inside OpenCode

API keys are read from the named environment variable, so they never sit in the config file.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

KINDS = ("anthropic", "openai", "opencode")
# OpenCode's gateway: every Zen text model is served through Chat Completions (translated where the
# model's native API differs); the Go subscription lives under /go and needs a per-chat session id.
OPENCODE = {
    "zen": dict(base_url="https://opencode.ai/inference/openai/v1"),
    "go": dict(
        base_url="https://opencode.ai/inference/go/openai/v1", session_header="x-opencode-session"
    ),
}


@dataclass
class Provider:
    name: str
    kind: str
    models: list[str]
    base_url: str = ""  # empty = the vendor's own API
    api_key_env: str = ""
    local: bool = False  # data stays on site: allowed for local-only projects
    vision: bool = True  # can look at rendered figures
    effort: str = ""  # anthropic: output_config.effort (low … max); empty = model default
    max_tokens: int = 16000
    stream_usage: bool = True  # openai: ask for token usage in the stream (stream_options)
    auth: str = ""  # "bearer": send the key as Authorization: Bearer (gateways that want it)
    headers: dict = field(default_factory=dict)  # extra request headers for a gateway
    session_header: str = ""  # header that carries the conversation id (e.g. OpenCode Go)
    plan: str = ""  # opencode: "zen" (default) or "go"
    # TLS trust: "" = Python's bundle; "system" = the OS store (networks that inspect HTTPS, e.g.
    # Cloudflare Gateway, install their root there); or the path of a PEM file
    ca_bundle: str = ""

    @property
    def official(self) -> bool:
        """The vendor's own API (vendor-only request features are sent only there)."""
        return not self.base_url

    def api_key(self) -> str:
        return os.environ.get(self.api_key_env, "") if self.api_key_env else ""

    def public(self) -> dict:
        return dict(name=self.name, models=self.models, local=self.local, vision=self.vision)


@dataclass
class AgentSettings:
    providers: list[Provider] = field(default_factory=list)
    max_steps: int = 12
    daily_tokens_per_user: int = 0
    local_only_projects: list[str] = field(default_factory=list)
    local_only_all: bool = False

    @property
    def enabled(self) -> bool:
        return bool(self.providers)

    def available(self, project_id: str) -> list[Provider]:
        """Providers this project may use: only local ones for local-only projects, and only those
        whose API key is present (a provider without api_key_env needs none)."""
        local_only = self.local_only_all or project_id in self.local_only_projects
        return [
            p
            for p in self.providers
            if (p.local or not local_only) and (not p.api_key_env or p.api_key())
        ]

    def missing_keys(self) -> list[str]:
        return [
            f"{p.name} ({p.api_key_env})"
            for p in self.providers
            if p.api_key_env and not p.api_key()
        ]

    def provider(self, project_id: str, name: str, model: str) -> Provider:
        for p in self.available(project_id):
            if p.name == name:
                if model not in p.models:
                    raise ValueError(f"{name} does not offer the model {model!r}")
                return p
        raise ValueError(f"The assistant {name!r} is not available for this project")


def agent_settings(raw: dict | None) -> AgentSettings:
    raw = dict(raw or {})
    providers = []
    for item in raw.pop("providers", None) or []:
        item = dict(item)
        if item.get("kind") not in KINDS:
            raise ValueError(f"agent provider {item.get('name')!r}: kind must be one of {KINDS}")
        if not item.get("name") or not item.get("models"):
            raise ValueError("every agent provider needs a name and a list of models")
        unknown = set(item) - set(Provider.__dataclass_fields__)
        if unknown:
            raise ValueError(f"agent provider {item['name']!r}: unknown keys {sorted(unknown)}")
        if item["kind"] == "opencode":
            item = _opencode(item)
        providers.append(Provider(**item))
    if len({p.name for p in providers}) != len(providers):
        raise ValueError("agent provider names must be unique")
    unknown = set(raw) - set(AgentSettings.__dataclass_fields__)
    if unknown:
        raise ValueError(f"unknown [agent] settings {sorted(unknown)}")
    return AgentSettings(providers=providers, **raw)


def _opencode(item):
    """Fill in OpenCode Zen / Go endpoints and checks; explicit settings in the config still win."""
    plan = item.setdefault("plan", "zen")
    if plan not in OPENCODE:
        raise ValueError(f"agent provider {item['name']!r}: plan must be 'zen' or 'go'")
    free = [m for m in item["models"] if str(m).endswith("-free")]
    if free:
        raise ValueError(
            f"agent provider {item['name']!r}: {', '.join(free)} only work inside OpenCode itself; "
            "use paid Zen models or the Go plan"
        )
    item.setdefault("api_key_env", "OPENCODE_API_KEY")
    for key, value in OPENCODE[plan].items():
        item.setdefault(key, value)
    return item
