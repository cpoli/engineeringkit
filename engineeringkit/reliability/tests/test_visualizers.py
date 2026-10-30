"""Smoke tests for engineeringkit.reliability.visualizers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit.reliability import weibull_hazard, weibull_reliability
from engineeringkit.reliability.visualizers import plot_hazard, plot_survival


def test_plot_survival_step_and_line():
    t = np.linspace(0, 5, 11)
    ax = plot_survival(t, weibull_reliability(t, 2, 1.5), label="Weibull", step=False)
    assert isinstance(ax, matplotlib.axes.Axes)
    plot_survival([0, 1, 2], [1, 0.5, 0.25], ax=ax)
    assert len(ax.lines) == 2
    with pytest.raises(ValueError):
        plot_survival([0, 1], [1.0])
    plt.close("all")


def test_plot_hazard_linear_and_log():
    t = np.linspace(0.1, 5, 11)
    ax = plot_hazard(t, weibull_hazard(t, 2, 0.5), label="beta 0.5", log=True)
    assert ax.get_yscale() == "log"
    assert plot_hazard(t, weibull_hazard(t, 2, 2)).get_yscale() == "linear"
    with pytest.raises(ValueError):
        plot_hazard([], [])
    plt.close("all")


def test_plot_weibull_probability_points_line_and_errors():
    from engineeringkit.reliability.visualizers import plot_weibull_probability

    ax = plot_weibull_probability([10, 20, 40], [0.2, 0.5, 0.8])
    assert len(ax.lines) == 1 and ax.get_legend() is None
    plot_weibull_probability([10, 20, 40], [0.2, 0.5, 0.8], eta=25, beta=1.5, ax=ax)
    x, y = ax.lines[-1].get_data()
    np.testing.assert_allclose(y, 1.5 * (x - np.log(25)))
    for times, f, kwargs in (([0, 1], [0.2, 0.5], {}), ([1, 2], [0.2, 1.0], {}), ([1, 2], [0.2, 0.5], {"eta": 3})):
        with pytest.raises(ValueError):
            plot_weibull_probability(times, f, **kwargs)
    plt.close("all")
