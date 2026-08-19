"""Small local signal benchmark."""

from pathlib import Path
from timeit import timeit

from genomecanvas.io.bedgraph import read_bedgraph


def main() -> None:
    path = Path(__file__).parents[1] / "tests/fixtures/mini_genome/signal_a.bedgraph"
    print(timeit(lambda: read_bedgraph(path, chrom="chr1", start=0, end=300), number=1000))


if __name__ == "__main__":
    main()
