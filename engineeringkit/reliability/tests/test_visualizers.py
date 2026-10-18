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
