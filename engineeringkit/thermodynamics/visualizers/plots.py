"""Plotting helpers for engineeringkit.thermodynamics: p-v and T-s diagrams."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics.systems.steam import saturation_pressure, steam_properties
from engineeringkit.thermodynamics.utils import if97

__all__ = ["plot_pv_diagram", "plot_ts_diagram", "plot_saturation_dome"]


def _xy(x, y, names):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or not x.size:
        raise ValueError(f"{names} must be matching nonempty vectors")
    return x, y


def plot_pv_diagram(v, p, ax=None, label=None):
    """Plot a process path on pressure-volume axes.

    Parameters
    ----------
    v : array_like
        Specific (or total) volumes along the path [m^3/kg or m^3].
    p : array_like
        Matching pressures [Pa], drawn in kPa.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    v, p = _xy(v, p, "v and p")
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(v, p / 1e3, label=label)
    ax.set(xlabel="volume", ylabel="pressure [kPa]")
    if label:
        ax.legend()
    return ax


def plot_saturation_dome(ax=None, T_min=273.16, T_max=647.096, n=120):
    """Draw the saturated-liquid and saturated-vapour lines of water on T-s axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    T_min, T_max : float, optional
        Temperature range [K] of the dome, within 273.16 K and the critical
        647.096 K.
    n : int, optional
        Number of points on each line.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if not 273.16 <= T_min < T_max <= 647.096:
        raise ValueError("require 273.16 <= T_min < T_max <= 647.096")
    if ax is None:
        _, ax = plt.subplots()
    # the region 4 equation overshoots 22.064 MPa at 647.096 K by a few parts in 1e11
    p_sat = np.minimum(saturation_pressure(np.linspace(T_min, T_max, n)), if97.P_CRIT)
    for x in (0.0, 1.0):
        states = [steam_properties(float(p), x=x) for p in p_sat]
        ax.plot([st.s / 1e3 for st in states], [st.T - 273.15 for st in states], "k-", lw=1)
    ax.set(xlabel="entropy [kJ/(kg K)]", ylabel="temperature [°C]")
    return ax


def plot_ts_diagram(s, T, ax=None, label=None, dome=True):
    """Plot a process path on temperature-entropy axes, optionally over the water dome.

    Parameters
    ----------
    s : array_like
        Specific entropies along the path [J/(kg K)], drawn in kJ/(kg K).
    T : array_like
        Matching temperatures [K], drawn in °C.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.
    dome : bool, optional
        Draw the saturation dome of water first.

    Returns
    -------
    matplotlib.axes.Axes
    """
    s, T = _xy(s, T, "s and T")
    if ax is None:
        _, ax = plt.subplots()
    if dome:
        plot_saturation_dome(ax)
    ax.plot(s / 1e3, T - 273.15, label=label)
    ax.set(xlabel="entropy [kJ/(kg K)]", ylabel="temperature [°C]")
    if label:
        ax.legend()
    return ax
