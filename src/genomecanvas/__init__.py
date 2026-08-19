"""GenomeCanvas: declarative genomic figure infrastructure."""

__version__ = "0.1.0"

from .model.background_spec import BackgroundSpec, TransformationSpec
from .model.figure_spec import FigureSpec
from .model.region import GenomeRegion
from .render.renderer import render
from .tracks.derived_signal import DerivedTrack

__all__ = [
    "BackgroundSpec", "DerivedTrack", "FigureSpec", "GenomeRegion",
    "TransformationSpec", "__version__", "render",
]
