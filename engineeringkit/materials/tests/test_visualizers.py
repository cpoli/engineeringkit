"""Smoke tests for engineeringkit.materials.visualizers and the result plots."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit import materials as m
from engineeringkit.materials.visualizers import plot_crack_growth, plot_sn_curve, plot_stress_strain, plot_yield_surfaces


def test_plot_stress_strain_scales_to_mpa():
    ax = plot_stress_strain([0.0, 0.001], [0.0, 200e6], label="steel")
    assert isinstance(ax, matplotlib.axes.Axes)
    assert ax.lines[0].get_ydata()[-1] == 200
    with pytest.raises(ValueError):
        plot_stress_strain([0.0, 0.001], [0.0])
    plt.close("all")


def test_plot_yield_surfaces_draws_the_hexagon_and_its_circumscribed_ellipse():
    ax = plot_yield_surfaces(250e6)
    tresca, mises = ax.lines[0], ax.lines[1]
    assert len(tresca.get_xdata()) == 7
    s1, s2 = mises.get_xdata(), mises.get_ydata()
    np.testing.assert_allclose(s1 * s1 - s1 * s2 + s2 * s2, 1)
    with pytest.raises(ValueError):
        plot_yield_surfaces(0.0)
    plt.close("all")


def test_plot_sn_curve_and_crack_growth():
    ax = plot_sn_curve([1e3, 1e6], [500e6, 250e6], label="fit")
    plot_sn_curve([1e4], [400e6], ax=ax, style="o")
    assert ax.get_xscale() == "log" and len(ax.lines) == 2
    with pytest.raises(ValueError):
        plot_sn_curve([0.0, 1e3], [1e6, 1e6])
    ax = plot_crack_growth([0.0, 100.0], [1e-3, 2e-3], label="a")
    assert ax.lines[0].get_ydata()[-1] == 2
    plt.close("all")


def test_result_plots_return_figure_and_axes():
    curve = m.StressStrainCurve([0.0, 0.001, 0.002], [0.0, 200e6, 250e6])
    fig, ax = curve.plot()
    assert ax.figure is fig
    growth = m.crack_growth_life(1e-29, 3.0, 100e6, 1e-3, 1e-2)
    fig, ax = growth.plot()
    assert ax.get_ylabel() == "crack length [mm]"
    plt.close("all")
