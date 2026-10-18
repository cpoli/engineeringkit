"""Smoke tests for engineeringkit.signal_processing.visualizers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pytest

from engineeringkit.signal_processing import butterworth_lowpass
from engineeringkit.signal_processing.visualizers import plot_magnitude_response, plot_pole_zero


def test_plot_magnitude_response_in_hz_and_rad_per_sample():
    b, a = butterworth_lowpass(4, 100.0, 1000.0)
    ax = plot_magnitude_response(b, a, fs=1000.0, label="order 4")
    assert isinstance(ax, matplotlib.axes.Axes)
    assert ax.get_xlabel() == "frequency [Hz]"
    sos = butterworth_lowpass(8, 100.0, 1000.0, output="sos")
    plot_magnitude_response(sos, None, fs=1000.0, ax=ax)
    assert len(ax.lines) == 2
    assert plot_magnitude_response([0.5, 0.5]).get_xlabel() == "frequency [rad/sample]"
    with pytest.raises(ValueError):
        plot_magnitude_response(b, a, fs=0.0)
    plt.close("all")


def test_plot_pole_zero_places_roots_in_z():
    ax = plot_pole_zero([1.0, 0.0, -1.0], [1.0, -0.5])
    zeros, poles = ax.lines[1], ax.lines[2]
    np.testing.assert_allclose(np.sort(zeros.get_xdata()), [-1.0, 1.0])
    np.testing.assert_allclose(np.sort(poles.get_xdata()), [0.0, 0.5])
    with pytest.raises(ValueError):
        plot_pole_zero([0.0], [1.0])
    plt.close("all")
