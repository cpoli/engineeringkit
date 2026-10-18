"""Plotting helpers for engineeringkit.signal_processing: magnitude responses and pole-zero maps."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit._validation import positive as _positive
from engineeringkit.signal_processing.systems.filters import frequency_response

__all__ = ["plot_magnitude_response", "plot_pole_zero"]


def plot_magnitude_response(b, a=1, fs=2 * np.pi, n=2048, ax=None, label=None):
    """Plot a digital filter's gain in decibels from DC to the Nyquist frequency.

    Parameters
    ----------
    b, a : array_like
        Filter polynomials, or second-order sections in ``b`` with ``a=None``,
        as in :func:`~engineeringkit.signal_processing.frequency_response`.
    fs : float, optional
        Sampling rate [Hz]; the default :math:`2\\pi` gives rad/sample.
    n : int, optional
        Number of frequencies.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    _positive(fs=fs)
    f, H = frequency_response(b, a, worN=n, fs=fs)
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(f, 20 * np.log10(np.maximum(np.abs(H), 1e-12)), label=label)
    ax.set(xlabel="frequency [Hz]" if fs != 2 * np.pi else "frequency [rad/sample]", ylabel="gain [dB]")
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()
    return ax


def plot_pole_zero(b, a, ax=None):
    """Plot the zeros (o) and poles (x) of :math:`B(z)/A(z)` with the unit circle.

    Parameters
    ----------
    b, a : array_like
        Numerator and denominator coefficients in :math:`z^{-1}` [-].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    b, a = np.atleast_1d(np.asarray(b, dtype=float)), np.atleast_1d(np.asarray(a, dtype=float))
    if b.ndim != 1 or a.ndim != 1 or not np.any(a) or not np.any(b):
        raise ValueError("b and a must be nonzero 1D coefficient vectors")
    # Pad to equal length so roots in z^-1 become roots in z, including those at the origin.
    n = max(b.size, a.size)
    zeros, poles = np.roots(np.pad(b, (0, n - b.size))), np.roots(np.pad(a, (0, n - a.size)))
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    theta = np.linspace(0, 2 * np.pi, 361)
    ax.plot(np.cos(theta), np.sin(theta), color="gray", lw=0.8)
    ax.plot(zeros.real, zeros.imag, "o", mfc="none", label="zeros")
    ax.plot(poles.real, poles.imag, "x", label="poles")
    ax.set_aspect("equal")
    ax.set(xlabel="real", ylabel="imaginary")
    ax.legend()
    return ax
