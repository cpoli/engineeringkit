"""Smoke tests for engineeringkit.communications.visualizers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np
import pytest

from engineeringkit.communications import ber_bpsk_awgn, qpsk
from engineeringkit.communications.visualizers import plot_ber_curve, plot_constellation


def test_plot_constellation_returns_axes():
    ax = plot_constellation(qpsk([0, 0, 0, 1, 1, 0, 1, 1]), label="QPSK")
    assert isinstance(ax, matplotlib.axes.Axes)
    with pytest.raises(ValueError):
        plot_constellation([])


def test_plot_ber_curve_returns_axes():
    x = np.linspace(0, 10, 11)
    ax = plot_ber_curve(x, ber_bpsk_awgn(10 ** (x / 10)), label="BPSK")
    assert isinstance(ax, matplotlib.axes.Axes)
    with pytest.raises(ValueError):
        plot_ber_curve([0, 1], [0.1])
