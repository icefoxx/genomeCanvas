"""Minimal CLI entry point for the scaffold release."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..model.figure_spec import FigureSpec
from ..render.renderer import render


def main() -> int:
    """Run the scaffold command-line interface."""
    parser = argparse.ArgumentParser(prog="gcanvas")
    parser.add_argument("command", nargs="?", choices=["render", "validate", "inspect", "init"])
    parser.add_argument("figure", nargs="?")
    parser.add_argument("--output", default="figure.pdf")
    arguments = parser.parse_args()
    if arguments.command == "validate":
        FigureSpec.from_yaml(arguments.figure)
        print("valid")
    elif arguments.command == "render":
        render(FigureSpec.from_yaml(arguments.figure), arguments.output)
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
