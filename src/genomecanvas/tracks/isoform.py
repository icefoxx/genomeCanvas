"""Individual transcript structure track."""

from .gene import GeneTrack


class IsoformTrack(GeneTrack):
    """Render individual isoforms using the same normalized transcript data."""

    type = "isoform"
