"""Plotting helpers for engineeringkit.fluid_engineering: the Moody chart and grade lines."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import vector as _vector

__all__ = ["plot_moody_diagram", "plot_grade_lines"]

DEFAULT_ROUGHNESS = (0.0, 1e-5, 1e-4, 1e-3, 5e-3, 0.01, 0.02, 0.05)


def plot_moody_diagram(relative_roughness=DEFAULT_ROUGHNESS, ax=None):
    """Draw a Moody chart: Darcy friction factor against Reynolds number.

    The laminar line :math:`64/Re` is drawn up to 2300 and one Colebrook
    curve per relative roughness from 4000 to :math:`10^8`; the transition
    range between them is shaded.

    Parameters
    ----------
    relative_roughness : array_like, optional
        Values of :math:`\\varepsilon/D` [-] in [0, 0.05], one curve each.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    from engineeringkit.fluid_engineering.systems.pipe_flow import colebrook_friction_factor, laminar_friction_factor

    roughness = _vector("relative_roughness", np.atleast_1d(relative_roughness))
    _interval("relative_roughness", roughness, 0, 0.05)
    if ax is None:
        _, ax = plt.subplots()
    Re_lam = np.geomspace(600, 2299, 30)
    ax.loglog(Re_lam, laminar_friction_factor(Re_lam), color="k", label="laminar, 64/Re")
    Re_turb = np.geomspace(4000, 1e8, 120)
    for eps in roughness:
        ax.loglog(Re_turb, colebrook_friction_factor(Re_turb, eps), label=f"ε/D = {eps:g}")
    ax.axvspan(2300, 4000, color="0.85", zorder=0)
    ax.set(xlabel="Reynolds number", ylabel="Darcy friction factor", title="Moody diagram")
    ax.grid(True, which="both", lw=0.3)
    ax.legend(fontsize="small")
    return ax


def plot_grade_lines(x, z, pressure_head, velocity_head, ax=None):
    """Draw the pipe centreline with its hydraulic and energy grade lines.

    The hydraulic grade line (HGL) is :math:`z + p/(\\rho g)`, the height
    water would rise in a piezometer; the energy grade line (EGL) adds the
    velocity head. Friction makes both fall in the direction of flow.

    Parameters
    ----------
    x : array_like
        Distance along the pipe [m], shape ``(n,)``.
    z : array_like
        Centreline elevation [m], shape ``(n,)``.
    pressure_head : array_like
        :math:`p/(\\rho g)` [m], shape ``(n,)``.
    velocity_head : array_like
        :math:`v^2/(2g)` [m], shape ``(n,)``.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x = _vector("x", x)
    z, hp, hv = (_vector(name, v) for name, v in (("z", z), ("pressure_head", pressure_head), ("velocity_head", velocity_head)))
    if not z.size == hp.size == hv.size == x.size:
        raise ValueError("x, z, pressure_head and velocity_head must have the same length")
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, z, color="0.4", lw=3, label="pipe centreline")
    ax.plot(x, z + hp, label="hydraulic grade line")
    ax.plot(x, z + hp + hv, "--", label="energy grade line")
    ax.set(xlabel="distance along pipe [m]", ylabel="elevation [m]")
    ax.legend()
    return ax
