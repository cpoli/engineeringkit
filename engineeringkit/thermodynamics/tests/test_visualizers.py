"""Smoke tests for engineeringkit.thermodynamics.visualizers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit.thermodynamics import isentropic_pressure_volume, steam_properties
from engineeringkit.thermodynamics.visualizers import plot_pv_diagram, plot_saturation_dome, plot_ts_diagram


def test_plot_pv_diagram_returns_axes():
    v = np.linspace(1.0, 0.125, 20)
    ax = plot_pv_diagram(v, isentropic_pressure_volume(1e5, 1.0, v), label="compression")
    assert isinstance(ax, matplotlib.axes.Axes)
    with pytest.raises(ValueError):
        plot_pv_diagram([1, 2], [1])
    plt.close("all")


def test_plot_ts_diagram_draws_the_dome_and_path():
    states = [steam_properties(10e3, x=0.0), steam_properties(10e6, 773.15), steam_properties(10e3, x=0.79)]
    ax = plot_ts_diagram([st.s for st in states], [st.T for st in states], label="Rankine")
    assert isinstance(ax, matplotlib.axes.Axes) and len(ax.lines) == 3
    assert len(plot_ts_diagram([0.0, 1.0], [300.0, 300.0], dome=False).lines) == 1
    _, ax = plt.subplots()
    assert plot_saturation_dome(ax, T_min=300.0, T_max=500.0, n=10) is ax
    liquid, vapour = plot_saturation_dome(n=5).lines
    # both lines end at the critical point, 373.946 °C
    assert np.allclose(liquid.get_xydata()[-1], vapour.get_xydata()[-1], rtol=1e-3)
    assert np.isclose(liquid.get_ydata()[-1], 373.946)
    with pytest.raises(ValueError):
        plot_saturation_dome(T_min=200.0)
    with pytest.raises(ValueError):
        plot_ts_diagram([], [])
    plt.close("all")
