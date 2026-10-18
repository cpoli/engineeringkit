"""Smoke tests for engineeringkit.control_systems.visualizers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit.control_systems import step_response, transfer_function
from engineeringkit.control_systems.visualizers import plot_bode, plot_nyquist, plot_response, plot_root_locus

G = transfer_function([1], [1, 2, 1])


def test_plot_response_returns_axes():
    t, y = step_response(G, np.linspace(0, 5, 51))
    ax = plot_response(t, y, label="step")
    assert isinstance(ax, matplotlib.axes.Axes) and ax.get_legend() is not None
    with pytest.raises(ValueError):
        plot_response([0, 1], [0])
    plt.close("all")


def test_plot_bode_returns_two_axes():
    ax_mag, ax_phase = plot_bode(G, np.logspace(-1, 1, 20), label="G")
    assert ax_mag.get_xscale() == "log" and ax_phase.get_legend() is None
    assert len(ax_mag.lines) == 1
    plt.close("all")


def test_plot_nyquist_draws_both_frequency_signs():
    ax = plot_nyquist(G, label="G")
    assert isinstance(ax, matplotlib.axes.Axes) and ax.get_legend() is not None
    plt.close("all")


def test_plot_root_locus_marks_poles_and_zeros():
    ax = plot_root_locus([1, 3], [1, 3, 2, 0], np.linspace(0, 10, 11))
    labels = [line.get_label() for line in ax.lines]
    assert "open-loop poles" in labels and "open-loop zeros" in labels
    ax = plot_root_locus([1], [1, 3, 2, 0], np.linspace(0, 10, 11))
    assert "open-loop zeros" not in [line.get_label() for line in ax.lines]
    plt.close("all")
