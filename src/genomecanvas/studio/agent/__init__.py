"""Studio chat agent: the GenomeCanvas skills behind a chat window, on Claude, OpenAI or local
models."""

from .loop import ChatAgent, ClientGone
from .providers import ProviderError
from .settings import AgentSettings, Provider, agent_settings

__all__ = [
    "AgentSettings",
    "ChatAgent",
    "ClientGone",
    "Provider",
    "ProviderError",
    "agent_settings",
]
