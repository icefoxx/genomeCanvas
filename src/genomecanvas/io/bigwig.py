"""Optional low-level BigWig adapter."""


def read_bigwig(path: str, chrom: str, start: int, end: int) -> object:
    """Query a BigWig region when pyBigWig is installed."""
    try:
        import pyBigWig  # type: ignore[import-not-found]
    except ImportError as error:
        raise RuntimeError("BigWig support requires the optional pyBigWig dependency") from error
    with pyBigWig.open(path) as handle:
        return handle.values(chrom, start, end, numpy=True)
