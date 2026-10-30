"""Smoke tests for engineeringkit.mechanical_engineering.visualizers and the result plots."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit import mechanical_engineering as m
from engineeringkit.mechanical_engineering.visualizers import plot_frequency_response, plot_linkage, plot_mode_shapes


def test_plot_frequency_response_hides_nonpositive_values():
    ax = plot_frequency_response([0.5, 1.0, 1.5], [1.3, np.inf, 0.8], label="test")
    assert isinstance(ax, matplotlib.axes.Axes)
    assert np.isnan(ax.lines[0].get_ydata()[1])
    with pytest.raises(ValueError):
        plot_frequency_response([0, 1], [1])
    plt.close("all")


def test_plot_mode_shapes_and_linkage():
    ax = plot_mode_shapes([[1.0, -0.5], [0.2, 1.0]], [10.0, 20.0])
    assert len(ax.lines) == 3
    with pytest.raises(ValueError):
        plot_mode_shapes([[1.0, -0.5]], [10.0, 20.0])
    ax = plot_linkage([[0, 0], [1, 1], [2, 0]], color="C1")
    assert isinstance(ax, matplotlib.axes.Axes)
    with pytest.raises(ValueError):
        plot_linkage([[0, 0]])
    plt.close("all")


def test_result_plots_return_figure_and_axes():
    fig, ax = m.tuned_mass_damper(10.0, 1e4, 0.1).plot()
    assert len(ax.lines) == 2
    fig, ax = m.jeffcott_response(10.0, 1e6, 1e-4, np.linspace(0, 600, 50), 100.0).plot()
    assert ax.get_xlabel() == "speed / critical speed"
    fig, ax = m.holzer_frequencies([1.0, 1.0, 1.0], [1.0, 1.0], 5.0).plot()
    assert len(ax.get_legend().get_texts()) == 2
    plt.close("all")


def test_curriculum_result_plots():
    M, K = np.diag([1.0, 2.0]), np.array([[300.0, -100.0], [-100.0, 100.0]])
    fig, ax = m.normal_modes(M, K).plot()
    assert ax.get_xlabel() == "degree of freedom"
    fig, ax = m.modal_response(M, K, [0.01, 0.0], t=np.linspace(0, 1, 11)).plot()
    assert len(ax.lines) == 2
    for load in (8e3, 1e5):  # in contact, and separated
        fig, ax = m.bolted_joint(1e9, 3e9, 20e3, load).plot()
        assert len(ax.lines) == 4
    fig, ax = m.balance_rotating_masses([1.0, 2.0], [0.1, 0.1], [0.0, 2.0], [0.1, 0.3], [0.0, 0.4]).plot()
    assert len(ax.get_legend().get_texts()) == 2
    fig, ax = m.fourbar_kinematics(1.0, 3.0, 2.0, 3.0, np.linspace(0, 2 * np.pi, 13), 1.0).plot()
    assert len(ax.lines) == 2
    fig, ax = m.slider_crank_kinematics(0.05, 0.2, np.linspace(0, 2 * np.pi, 13), 100.0).plot()
    assert ax.get_ylabel() == "slider velocity [m/s]"
    fig, axs = m.cam_follower_motion(np.linspace(0, 1, 11), 0.01, 1.0).plot()
    assert len(axs) == 4 and axs[0].get_title() == "cycloidal rise"
    plt.close("all")
