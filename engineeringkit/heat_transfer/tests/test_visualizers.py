"""Smoke tests for engineeringkit.heat_transfer.visualizers and result plots."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit.heat_transfer import conduction_1d
from engineeringkit.heat_transfer.visualizers import plot_blackbody_spectra, plot_temperature_profiles


def test_plot_temperature_profiles():
    x = np.linspace(0, 1, 5)
    ax = plot_temperature_profiles(x, [300, 310, 320, 310, 300])
    assert isinstance(ax, matplotlib.axes.Axes)
    ax = plot_temperature_profiles(x, np.full((2, 5), 300.0), times=[0, 1])
    assert len(ax.get_legend().get_texts()) == 2
    with pytest.raises(ValueError):
        plot_temperature_profiles(x, [300, 310])
    with pytest.raises(ValueError):
        plot_temperature_profiles(x, np.full((2, 5), 300.0), times=[0])
    plt.close("all")


def test_plot_blackbody_spectra():
    ax = plot_blackbody_spectra([1000.0, 5772.0])
    assert len(ax.get_lines()) == 4
    ax = plot_blackbody_spectra(300.0, wavelength=np.geomspace(1e-6, 1e-4, 50), ax=ax)
    with pytest.raises(ValueError):
        plot_blackbody_spectra([0.0])
    plt.close("all")


def test_conduction_result_plot():
    result = conduction_1d([300.0] * 6, 0.1, 1e-5, 1.0, 10, T_left=400.0)
    fig, ax = result.plot(profiles=3)
    assert len(ax.get_lines()) == 3
    plt.close(fig)


def test_conduction_2d_and_enclosure_results():
    from engineeringkit.heat_transfer import conduction_2d, radiation_enclosure

    result = conduction_2d(1.0, 0.5, 9, 5, 300.0, 310.0, 320.0, 400.0)
    summary = result.summary()
    assert summary["nodes"] == 45 and summary["max_K"] == 400.0
    assert abs(summary["imbalance_W_per_m"]) < 1e-9
    fig, ax = result.plot(levels=5)
    assert ax.get_xlabel() == "x [m]"
    plt.close(fig)
    enclosure = radiation_enclosure([500.0, 300.0], [0.8, 0.8], [1.0, 1.0], [[0.0, 1.0], [1.0, 0.0]])
    assert abs(enclosure.summary()["imbalance_W"]) < 1e-9
    assert enclosure.summary()["temperature_K"] == [500.0, 300.0]
    fig, ax = enclosure.plot()
    assert len(ax.patches) == 2
    plt.close(fig)
