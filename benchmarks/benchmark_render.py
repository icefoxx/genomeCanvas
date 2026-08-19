"""Small local rendering benchmark."""

from pathlib import Path
from timeit import timeit

from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.render.renderer import render


def main() -> None:
    root = Path(__file__).parents[1]
    spec = FigureSpec.from_yaml(root / "examples/basic_signal/figure.yaml")
    print(timeit(lambda: render(spec, "/var/tmp/genomecanvas-benchmark.pdf"), number=1))


if __name__ == "__main__":
    main()
