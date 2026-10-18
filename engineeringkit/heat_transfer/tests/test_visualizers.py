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
