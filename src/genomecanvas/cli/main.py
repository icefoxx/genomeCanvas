"""Minimal CLI entry point for the scaffold release."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..model.figure_spec import FigureSpec
from ..render.renderer import render


def main(argv: list[str] | None = None) -> int:
    """Run the command-line interface."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "heatmap":
        from .heatmap import run_heatmap

        return run_heatmap(argv[1:])
    parser = argparse.ArgumentParser(prog="gcanvas")
    parser.add_argument(
        "command", nargs="?", choices=["render", "validate", "inspect", "init", "heatmap"]
    )
    parser.add_argument("figure", nargs="?")
    parser.add_argument("--output", default="figure.pdf")
    parser.add_argument(
        "--bigwig-engine",
        choices=["auto", "pybigwig", "pybigtools", "bwtool"],
        default=None,
        help="BigWig backend (overrides figure YAML and GENOMECANVAS_BIGWIG_ENGINE)",
    )
    arguments = parser.parse_args(argv)
    if arguments.command == "validate":
        FigureSpec.from_yaml(arguments.figure)
        print("valid")
    elif arguments.command == "render":
        spec = FigureSpec.from_yaml(arguments.figure)
        if arguments.bigwig_engine:
            spec = spec.copy(update={"bigwig_engine": arguments.bigwig_engine})
        render(spec, arguments.output)
    elif arguments.command == "inspect":
        from .inspect import inspect_file

        print(inspect_file(arguments.figure))
    elif arguments.command == "init":
        target = Path(arguments.figure or "figure.yaml")
        target.write_text(
            "version: '1.0'\n"
            "genome:\n  assembly: hg38\n"
            "region:\n  chrom: chr1\n  start: 0\n  end: 1000\n"
            "tracks: []\n",
            encoding="utf-8",
        )
        print(target)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
