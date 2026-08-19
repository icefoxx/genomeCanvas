from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from genomecanvas.render.plotting import make_grid, plot_signal, save_figure


def test_plotting_helpers_create_and_export_standard_frame(tmp_path: Path) -> None:
    figure, axes = make_grid(2, 2, figsize=(4, 3), sharex="col")
    x = np.arange(4)
    plot_signal(axes[0, 0], x, x.astype(float), color="#3366aa")

    save_figure(figure, tmp_path, "example")

    assert axes.shape == (2, 2)
    assert (tmp_path / "example.png").is_file()
    assert (tmp_path / "example.pdf").is_file()
    assert not plt.fignum_exists(figure.number)
