"""Core genomic data structures."""

from .coordinates import CoordinateConversion, chromosome_alias
from .intervals import GenomicInterval

__all__ = ["CoordinateConversion", "GenomicInterval", "chromosome_alias"]
