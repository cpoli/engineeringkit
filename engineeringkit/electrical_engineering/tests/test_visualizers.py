"""Smoke tests for engineeringkit.electrical_engineering.visualizers and result plots."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit import electrical_engineering as e
from engineeringkit.electrical_engineering.visualizers import plot_phasors, plot_smith_chart


def test_plot_phasors_returns_axes():
    ax = plot_phasors(e.phase_components(0, 1, 0), labels=["a", "b", "c"])
    assert isinstance(ax, matplotlib.axes.Axes)
    assert isinstance(plot_phasors([0j]), matplotlib.axes.Axes)
    for bad in ([], [np.nan]):
        with pytest.raises(ValueError):
            plot_phasors(bad)
    with pytest.raises(ValueError):
        plot_phasors([1, 1j], labels=["a"])
    plt.close("all")


def test_plot_smith_chart_returns_axes():
    assert isinstance(plot_smith_chart(), matplotlib.axes.Axes)
    gamma = e.reflection_coefficient(e.input_impedance(30 + 40j, 50, np.linspace(0, np.pi, 50)), 50)
    assert isinstance(plot_smith_chart(gamma, label="line"), matplotlib.axes.Axes)
    assert isinstance(plot_smith_chart(0.5), matplotlib.axes.Axes)
    with pytest.raises(ValueError):
        plot_smith_chart([1.5])
    plt.close("all")


def test_result_plots():
    nodal = e.nodal_analysis(1, [(1, 0, 1.0)], current_sources=[(0, 1, 1.0)])
    flow = e.gauss_seidel_load_flow(e.admittance_matrix(2, [(0, 1, 0.1j)]), [0.5], [0.1])
    mesh = e.mesh_analysis(1, [(1, 0, 2.0)], voltage_sources=[(1, 0, 4.0)])
    swing = e.swing_response(5.0, 1.0, 2.0, 0.0, 2.0, 0.1, 1.0)
    assert swing.summary() == {"stable": True, "max_delta": float(swing.delta.max()), "clearing_time": 0.1}
    assert mesh.summary()["max_current"] == 2.0
    for result in (nodal, flow, mesh, swing, e.SeriesRLC(2, 1, 0.25).step_response()):
        fig, ax = result.plot()
        assert isinstance(ax, matplotlib.axes.Axes)
    _, ax = plt.subplots()
    assert nodal.plot(ax=ax)[1] is ax
    plt.close("all")
