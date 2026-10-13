"""Plotting helpers for engineeringkit.aerodynamics: sections, polars and the atmosphere."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics.core.base import Airfoil, AtmosphereState

__all__ = ["plot_airfoil", "plot_drag_polar", "plot_atmosphere"]


def plot_airfoil(airfoil: Airfoil, ax=None):
    """Draw an airfoil outline and its mean camber line on equal axes.

    Parameters
    ----------
    airfoil : Airfoil
        From :func:`engineeringkit.aerodynamics.naca4_airfoil`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 2.5))
    outline = np.vstack((airfoil.upper[::-1], airfoil.lower[1:]))
    ax.fill(outline[:, 0], outline[:, 1], alpha=0.25)
    ax.plot(outline[:, 0], outline[:, 1], "k", lw=1)
    ax.plot(airfoil.x, airfoil.camber, "C3--", lw=1, label="mean camber line")
    ax.set_aspect("equal")
    ax.set(xlabel="x / c", ylabel="y / c", title=airfoil.name)
    ax.legend(loc="upper right")
    return ax


def plot_drag_polar(cl, cd, ax=None, label=None):
    """Plot lift coefficient against drag coefficient and mark the best-glide tangent.

    Parameters
    ----------
    cl, cd : array_like
        Matching lift and drag coefficients [-].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label for the polar.

    Returns
    -------
    matplotlib.axes.Axes
    """
    cl, cd = np.asarray(cl, float), np.asarray(cd, float)
    if cl.shape != cd.shape or cl.ndim != 1 or not cl.size:
        raise ValueError("cl and cd must be matching nonempty vectors")
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(cd, cl, label=label)
    i = int(np.argmax(np.where(cd > 0, cl / np.where(cd > 0, cd, 1), -np.inf)))
    ax.plot([0, cd[i]], [0, cl[i]], ":", color=ax.lines[-1].get_color())
    ax.plot(cd[i], cl[i], "o", color=ax.lines[-1].get_color())
    ax.set(xlabel="drag coefficient $C_D$", ylabel="lift coefficient $C_L$")
    if label:
        ax.legend()
    return ax


def plot_atmosphere(state: AtmosphereState, axes=None):
    """Plot temperature, pressure and density against altitude side by side.

    Parameters
    ----------
    state : AtmosphereState
        From :func:`engineeringkit.aerodynamics.isa_atmosphere` with an altitude array.
    axes : sequence of three matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    ndarray of matplotlib.axes.Axes
    """
    if axes is None:
        _, axes = plt.subplots(1, 3, figsize=(10, 3.5), sharey=True)
    km = np.atleast_1d(state.altitude) / 1e3
    for ax, values, xlabel in zip(
        axes,
        (state.temperature, state.pressure / 1e3, state.density),
        ("temperature [K]", "pressure [kPa]", "density [kg/m$^3$]"),
        strict=True,
    ):
        ax.plot(np.atleast_1d(values), km)
        ax.set_xlabel(xlabel)
    axes[0].set_ylabel("altitude [km]")
    return np.asarray(axes)
