"""Track lifecycle protocol."""

from typing import Protocol


class TrackRenderer(Protocol):
    """Backend-independent lifecycle expected of renderable tracks."""

    def prepare(self) -> None: ...

    def measure(self) -> float: ...
