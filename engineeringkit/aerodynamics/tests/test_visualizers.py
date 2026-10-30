"""Smoke tests for engineeringkit.aerodynamics.visualizers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np
import pytest

from engineeringkit.aerodynamics import drag_polar, isa_atmosphere, naca4_airfoil
from engineeringkit.aerodynamics.visualizers import plot_airfoil, plot_atmosphere, plot_drag_polar


def test_plot_airfoil_returns_axes_and_result_plot_returns_figure():
    foil = naca4_airfoil("4412")
    assert isinstance(plot_airfoil(foil), matplotlib.axes.Axes)
    fig, ax = foil.plot()
    assert ax.figure is fig


def test_plot_drag_polar_marks_best_glide():
    cl = np.linspace(-0.5, 1.5, 81)
    ax = plot_drag_polar(cl, drag_polar(cl, 0.02, 8, 0.8), label="AR 8")
    assert isinstance(ax, matplotlib.axes.Axes)
    with pytest.raises(ValueError):
        plot_drag_polar([0.1, 0.2], [0.01])


def test_plot_atmosphere_returns_three_axes():
    state = isa_atmosphere(np.linspace(0, 20000, 41))
    axes = plot_atmosphere(state)
    assert len(axes) == 3
    fig, axes = state.plot()
    assert axes[0].figure is fig


def test_new_result_plots_and_summaries():
    from engineeringkit.aerodynamics import fanno_flow, lifting_line, oblique_shock, panel_method, rayleigh_flow
    from engineeringkit.aerodynamics.visualizers import plot_theta_beta_mach

    shock = oblique_shock(2.5, np.radians(12.0))
    fig, ax = shock.plot()
    assert ax.figure is fig and shock.summary()["mach2"] == shock.mach2
    assert isinstance(plot_theta_beta_mach([1.5, 3.0]), matplotlib.axes.Axes)
    M = np.linspace(0.3, 3.0, 50)
    for state in (fanno_flow(M), rayleigh_flow(M)):
        fig, ax = state.plot()
        assert ax.figure is fig
        assert state.summary()["mach"] == M[0]
    wing = lifting_line(8.0, 1.0, 0.1)
    fig, ax = wing.plot()
    assert ax.figure is fig and wing.summary()["CL"] == wing.CL
    panels = panel_method(naca4_airfoil("2412", 41), 0.05)
    fig, ax = panels.plot()
    assert ax.figure is fig and panels.summary()["panels"] == 80
