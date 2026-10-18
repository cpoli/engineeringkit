"""Plotting helpers for engineeringkit.heat_transfer: temperature profiles and blackbody spectra."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector

__all__ = ["plot_temperature_profiles", "plot_blackbody_spectra"]


def plot_temperature_profiles(x, T, times=None, ax=None):
    """Plot one or more temperature profiles :math:`T(x)`.

    Parameters
    ----------
    x : array_like
        Positions [m], shape ``(n,)``.
    T : array_like
        Temperatures [K], shape ``(n,)`` or ``(k, n)`` for ``k`` profiles.
    times : array_like, optional
        Time of each profile [s], shape ``(k,)``, used as legend labels.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x = _vector("x", x)
    T = np.atleast_2d(np.asarray(T, dtype=float))
    if T.ndim != 2 or T.shape[1] != x.size:
        raise ValueError("T must have shape (n,) or (k, n) with n = len(x)")
    if times is not None and np.size(times) != T.shape[0]:
        raise ValueError("times must give one time per profile")
    if ax is None:
        _, ax = plt.subplots()
    for i, row in enumerate(T):
        ax.plot(x, row, label=None if times is None else f"t = {np.ravel(times)[i]:.3g} s")
    ax.set(xlabel="position [m]", ylabel="temperature [K]")
    if times is not None:
        ax.legend()
    return ax


def plot_blackbody_spectra(temperatures, wavelength=None, ax=None):
    """Plot Planck's spectral emissive power for several temperatures, with each Wien peak marked.

    Parameters
    ----------
    temperatures : array_like
        Absolute temperatures [K], positive.
    wavelength : array_like, optional
        Wavelengths [m]; defaults to 0.1-100 um on a log grid.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    from engineeringkit.heat_transfer.systems.radiation import blackbody_peak_wavelength, planck_spectral_emissive_power

    temperatures = _vector("temperatures", np.atleast_1d(temperatures))
    _positive(temperatures=temperatures)
    lam = np.geomspace(1e-7, 1e-4, 600) if wavelength is None else _vector("wavelength", wavelength)
    if ax is None:
        _, ax = plt.subplots()
    for T in temperatures:
        line = ax.loglog(lam * 1e6, planck_spectral_emissive_power(lam, T), label=f"{T:g} K")[0]
        peak = blackbody_peak_wavelength(T)
        ax.plot(peak * 1e6, planck_spectral_emissive_power(peak, T), "o", color=line.get_color())
    ax.set(xlabel="wavelength [um]", ylabel="spectral emissive power [W/(m² m)]")
    ax.legend()
    return ax
