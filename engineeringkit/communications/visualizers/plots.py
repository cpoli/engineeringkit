"""Plotting helpers for engineeringkit.communications: constellations, error-rate curves, eye diagrams and antenna patterns."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_constellation", "plot_ber_curve", "plot_eye_diagram", "plot_antenna_pattern"]


def plot_constellation(symbols, ax=None, label=None):
    """Scatter complex symbols in the in-phase/quadrature plane.

    Parameters
    ----------
    symbols : array_like
        Complex (or real) symbols.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    s = np.asarray(symbols, dtype=complex).ravel()
    if not s.size:
        raise ValueError("symbols must be nonempty")
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    ax.plot(s.real, s.imag, ".", label=label)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_aspect("equal")
    ax.set(xlabel="in-phase", ylabel="quadrature")
    if label:
        ax.legend()
    return ax


def plot_ber_curve(ebn0_db, ber, ax=None, label=None):
    """Semilog plot of bit-error rate against :math:`E_b/N_0` in decibels.

    Parameters
    ----------
    ebn0_db : array_like
        :math:`E_b/N_0` values [dB].
    ber : array_like
        Matching bit-error rates [-]; zeros are omitted from the log axis.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, y = np.asarray(ebn0_db, dtype=float), np.asarray(ber, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or not x.size:
        raise ValueError("ebn0_db and ber must be matching nonempty vectors")
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogy(x, np.where(y > 0, y, np.nan), label=label)
    ax.set(xlabel="$E_b/N_0$ [dB]", ylabel="bit-error rate")
    ax.grid(True, which="both", alpha=0.3)
    if label:
        ax.legend()
    return ax


def plot_eye_diagram(signal, samples_per_symbol, ax=None, spans=2, max_traces=200, offset=0):
    """Overlay successive symbol-length slices of a waveform to show its eye opening.

    Parameters
    ----------
    signal : array_like
        Real baseband waveform sampled ``samples_per_symbol`` times per symbol.
    samples_per_symbol : int
        Oversampling factor, at least 2.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    spans : int, optional
        Symbol periods per trace, at least 1.
    max_traces : int, optional
        Largest number of traces drawn, at least 1.
    offset : int, optional
        Samples skipped before the first trace, for example a filter delay
        that puts a sampling instant at the centre of the plot.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x = np.asarray(signal, dtype=float)
    for name, value, minimum in (("samples_per_symbol", samples_per_symbol, 2), ("spans", spans, 1), ("max_traces", max_traces, 1)):
        if not isinstance(value, (int, np.integer)) or isinstance(value, bool) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    if not isinstance(offset, (int, np.integer)) or isinstance(offset, bool) or offset < 0:
        raise ValueError("offset must be a nonnegative integer")
    length = spans * samples_per_symbol + 1
    starts = np.arange(offset, x.size - length + 1, samples_per_symbol)[:max_traces]
    if x.ndim != 1 or not starts.size:
        raise ValueError("signal must be a 1-D waveform longer than one trace")
    if ax is None:
        _, ax = plt.subplots()
    t = np.arange(length) / samples_per_symbol
    ax.plot(t, np.stack([x[k : k + length] for k in starts]).T, color="C0", lw=0.6, alpha=0.4)
    ax.set(xlabel="time [symbol periods]", ylabel="amplitude", xlim=(0, spans))
    ax.grid(True, alpha=0.3)
    return ax


def plot_antenna_pattern(theta, pattern, ax=None, label=None, db_range=None):
    """Polar plot of an antenna's radiation pattern against the angle from its axis.

    Parameters
    ----------
    theta : array_like
        Angles [rad].
    pattern : array_like
        Matching nonnegative pattern values (field or power), normalised or not.
    ax : matplotlib.axes.Axes, optional
        Polar axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.
    db_range : float, optional
        If given, plot the power pattern in decibels relative to its peak,
        clipped this many dB below it.

    Returns
    -------
    matplotlib.axes.Axes
    """
    th, p = np.asarray(theta, dtype=float), np.asarray(pattern, dtype=float)
    if th.shape != p.shape or th.ndim != 1 or not th.size or np.any(p < 0) or not np.any(p > 0):
        raise ValueError("theta and pattern must be matching nonempty vectors of nonnegative values")
    if ax is None:
        _, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    if db_range is None:
        ax.plot(th, p, label=label)
    else:
        if not np.isfinite(db_range) or db_range <= 0:
            raise ValueError("db_range must be a positive number of decibels")
        level = 10 * np.log10(np.maximum(p / p.max(), 10 ** (-db_range / 10)))
        ax.plot(th, level + db_range, label=label)
        ticks = np.linspace(0, db_range, 4)
        ax.set_yticks(ticks, [f"{t - db_range:.0f} dB" for t in ticks])
    if label:
        ax.legend()
    return ax
