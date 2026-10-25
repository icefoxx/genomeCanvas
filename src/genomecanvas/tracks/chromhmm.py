"""Chromatin-state track backed by a 4-column chromHMM BED."""

from pathlib import Path

from ..io.chromhmm import ChromatinState, read_chromhmm


class ChromHMMTrack:
    """Load chromHMM intervals without interpreting biological state meaning."""

    type = "chromhmm"

    def __init__(self, source: str) -> None:
        self.source = source

    def states(self) -> list[ChromatinState]:
        """Return chromatin-state intervals from ``source``."""
        return read_chromhmm(Path(self.source))
