"""Region-level summaries, PCA, and simple differential peak tests."""

from .diffpeak import diffpeak
from .pca import pca_plot, pca_scores
from .summary import extract_summary
from .volcano import volcano_plot

__all__ = ["diffpeak", "extract_summary", "pca_plot", "pca_scores", "volcano_plot"]
