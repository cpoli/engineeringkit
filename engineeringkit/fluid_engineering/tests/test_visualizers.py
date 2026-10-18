"""Smoke tests for engineeringkit.fluid_engineering.visualizers and result plots."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit.fluid_engineering import buckingham_pi, hardy_cross
from engineeringkit.fluid_engineering.visualizers import plot_grade_lines, plot_moody_diagram


def test_plot_moody_diagram():
    ax = plot_moody_diagram()
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.get_lines()) == 9
    ax = plot_moody_diagram(0.001, ax=None)
    assert len(ax.get_lines()) == 2
    with pytest.raises(ValueError):
        plot_moody_diagram([0.1])
    plt.close("all")


def test_plot_grade_lines():
    x = np.linspace(0, 100, 5)
    ax = plot_grade_lines(x, np.zeros(5), np.linspace(10, 8, 5), np.full(5, 0.2))
    assert len(ax.get_lines()) == 3
    assert np.allclose(ax.get_lines()[2].get_ydata(), np.linspace(10.2, 8.2, 5))
    with pytest.raises(ValueError):
        plot_grade_lines(x, np.zeros(4), np.zeros(5), np.zeros(5))
    plt.close("all")


def test_result_plots():
    fig, ax = hardy_cross([1.0, 4.0], [[1, -1]], [0.15, 0.15]).plot()
    assert ax.get_yscale() == "log"
    plt.close(fig)
    pi = buckingham_pi({"rho": {"M": 1, "L": -3}, "V": {"L": 1, "T": -1}, "D": {"L": 1}, "mu": {"M": 1, "L": -1, "T": -1}})
    fig, ax = pi.plot()
    assert len(ax.texts) == 4
    plt.close(fig)
