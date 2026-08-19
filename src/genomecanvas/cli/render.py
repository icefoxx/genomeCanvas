"""Render command implementation."""

from ..model.figure_spec import FigureSpec
from ..render.renderer import render


def run(path: str, output: str = "figure.pdf") -> None:
    render(FigureSpec.from_yaml(path), output)
