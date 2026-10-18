import matplotlib.pyplot as plt
import pytest

from engineeringkit.visualizers import animate_displacement, plot_history


def test_scalar_history_plot():
    fig, ax = plot_history([0, 1], [1, 0])
    assert len(ax.lines) == 1
    plt.close(fig)


@pytest.mark.parametrize("times,values,labels", [([0, 1], [1], None), ([0, 1], [[1], [2]], ["x", "y"])])
def test_history_shape_and_labels(times, values, labels):
    with pytest.raises(ValueError):
        plot_history(times, values, labels)


@pytest.mark.parametrize("times,values,interval", [([0, 1], [1], 30), ([0, 1], [1, float("nan")], 30), ([0, 1], [0, 1], 0)])
def test_invalid_animation(times, values, interval):
    with pytest.raises(ValueError):
        animate_displacement(times, values, interval=interval)
