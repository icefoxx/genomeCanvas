"""Small local GTF benchmark."""

from pathlib import Path
from timeit import timeit

from genomecanvas.io.gtf import read_gtf


def main() -> None:
    path = Path(__file__).parents[1] / "tests/fixtures/mini_genome/mini.gtf"
    print(timeit(lambda: read_gtf(path), number=1000))


if __name__ == "__main__":
    main()
