"""Smoke tests for engineeringkit.thermodynamics.visualizers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit.thermodynamics import (
    brayton_cycle,
    isentropic_pressure_volume,
    rankine_cycle,
    steam_properties,
    vapor_compression_cycle,
)
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


def test_cycle_results_plot_their_state_diagrams():
    res = rankine_cycle(15e6, 10e3, 873.15, reheat_pressure=4e6, T_reheat=873.15, feedwater_pressure=1.2e6)
    fig, ax = res.plot()
    # dome (2 lines), one segment per process of the main loop, and the state markers
    assert len(ax.lines) == 2 + (len(res.labels) - 1) + 1 and ax.get_title().startswith("Rankine")
    gas = brayton_cycle(8.0, 300.0, 1300.0, regenerator_effectiveness=0.7, compressor_stages=2)
    _, ax = gas.plot()
    assert len(gas.labels) == 8 and len(ax.lines) == len(gas.labels) + 1  # a closed loop of 8 processes and the markers
    _, given = plt.subplots()
    fig, ax = vapor_compression_cycle(0.14e6, 0.8e6).plot(ax=given)
    assert ax is given and ax.get_yscale() == "log" and len(ax.lines) == 3
    assert len(vapor_compression_cycle(0.2e6, 1e6).plot()[1].lines) == 3
    plt.close("all")
