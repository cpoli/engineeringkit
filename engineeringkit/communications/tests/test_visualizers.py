"""Smoke tests for engineeringkit.communications.visualizers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np
import pytest

from engineeringkit.communications import ber_bpsk_awgn, bpsk, qpsk
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


def test_plot_eye_diagram_returns_axes():
    from engineeringkit.communications.visualizers import plot_eye_diagram

    signal = np.repeat(bpsk(np.random.default_rng(0).integers(0, 2, 40)), 8)
    ax = plot_eye_diagram(signal, 8, offset=4)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 38
    for kwargs in ({"samples_per_symbol": 1}, {"samples_per_symbol": 8, "spans": True}, {"samples_per_symbol": 8, "offset": -1}):
        with pytest.raises(ValueError):
            plot_eye_diagram(signal, **kwargs)
    with pytest.raises(ValueError):
        plot_eye_diagram(signal[:10], 8, spans=2)


def test_plot_antenna_pattern_linear_and_db():
    from engineeringkit.communications import half_wave_dipole_pattern
    from engineeringkit.communications.visualizers import plot_antenna_pattern

    theta = np.linspace(0, 2 * np.pi, 361)
    pattern = half_wave_dipole_pattern(theta) ** 2
    ax = plot_antenna_pattern(theta, pattern, label="dipole")
    assert ax.name == "polar"
    ax = plot_antenna_pattern(theta, pattern, db_range=30)
    assert min(ax.lines[0].get_ydata()) == pytest.approx(0)
    with pytest.raises(ValueError):
        plot_antenna_pattern(theta, pattern, db_range=0)
    with pytest.raises(ValueError):
        plot_antenna_pattern(theta, -pattern)


def test_huffman_code_plot():
    from engineeringkit.communications import huffman_code

    fig, ax = huffman_code([0.5, 0.3, 0.2, 0.0]).plot()
    assert ax.figure is fig and len(ax.patches) == 4 + 3
