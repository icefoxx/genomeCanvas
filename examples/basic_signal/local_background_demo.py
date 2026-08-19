"""Synthetic TASK-031 demonstration."""

import numpy as np

from genomecanvas.model.background_spec import BackgroundSpec
from genomecanvas.transform.background import estimate_background, transform_values


def main() -> None:
    background = np.full(100, 10.0)
    background[0] = 110.0
    estimate = estimate_background(
        background, background, BackgroundSpec(upstream=100, downstream=100)
    )
    peak = np.array([110.0])
    print("median background:", estimate.value)
    print("subtract:", transform_values(peak, estimate.value, "subtract"))
    print("ratio:", transform_values(peak, estimate.value, "ratio", 0.0))
    print("log2 ratio:", transform_values(peak, estimate.value, "log2_ratio", 0.0))


if __name__ == "__main__":
    main()
