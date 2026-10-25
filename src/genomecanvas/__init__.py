"""GenomeCanvas: declarative genomic figure infrastructure."""

__version__ = "0.1.0"

from .analysis import diffpeak, extract_summary, pca_plot, pca_scores, volcano_plot
from .interact import FigureEdit, OverlayAnnotation, apply_edit
from .io.bigwig import resolve_engine as resolve_bigwig_engine
from .io.bigwig import set_bigwig_engine, using_bigwig_engine
from .model.background_spec import BackgroundSpec, TransformationSpec
from .model.figure_spec import FigureSpec
from .model.region import GenomeRegion
from .pipeline import HeatmapPlotStyle, ProfilePlotStyle, run_heatmap_pipeline
from .profile import extract_profile, profile_heatmap, profile_plot, profile_summarize
from .render.renderer import render
from .tracks.derived_signal import DerivedTrack
from .transform.summarize import summarize_signals

__all__ = [
    "BackgroundSpec",
    "DerivedTrack",
    "FigureEdit",
    "FigureSpec",
    "GenomeRegion",
    "HeatmapPlotStyle",
    "ProfilePlotStyle",
    "TransformationSpec",
    "OverlayAnnotation",
    "__version__",
    "apply_edit",
    "diffpeak",
    "extract_profile",
    "extract_summary",
    "pca_plot",
    "pca_scores",
    "profile_heatmap",
    "profile_plot",
    "profile_summarize",
    "render",
    "run_heatmap_pipeline",
    "resolve_bigwig_engine",
    "set_bigwig_engine",
    "summarize_signals",
    "using_bigwig_engine",
    "volcano_plot",
]
