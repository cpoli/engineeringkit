"""Plotting helpers for engineeringkit.aerodynamics: sections, polars, atmosphere, shocks, duct flow, span loads, pressure."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics.core.base import Airfoil, AtmosphereState, DuctFlowState, LiftingLineResult, PanelResult

__all__ = [
    "plot_airfoil",
    "plot_drag_polar",
    "plot_atmosphere",
    "plot_theta_beta_mach",
    "plot_duct_flow",
    "plot_spanwise_loading",
    "plot_pressure_coefficient",
]


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


def plot_theta_beta_mach(mach_numbers, gamma=1.4, ax=None):
    """Draw the oblique-shock theta-beta-M chart: wave angle against deflection for each Mach number.

    Parameters
    ----------
    mach_numbers : sequence of float
        Upstream Mach numbers [-], each greater than one.
    gamma : float, optional
        Ratio of specific heats [-].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
        Weak branches are solid and strong branches dashed; the dotted line
        joins the maximum-deflection points.
    """
    from engineeringkit.aerodynamics.systems.supersonic import oblique_shock_deflection, oblique_shock_max_deflection

    if ax is None:
        _, ax = plt.subplots()
    peaks = []
    for M in mach_numbers:
        theta_max, beta_max = oblique_shock_max_deflection(M, gamma)
        beta = np.linspace(np.arcsin(1 / M), np.pi / 2, 200)
        theta = oblique_shock_deflection(M, beta, gamma)
        weak = beta <= beta_max
        (line,) = ax.plot(np.degrees(theta[weak]), np.degrees(beta[weak]), label=f"M = {M:g}")
        ax.plot(np.degrees(theta[~weak]), np.degrees(beta[~weak]), "--", color=line.get_color())
        peaks.append((np.degrees(theta_max), np.degrees(beta_max)))
    ax.plot(*np.array(peaks).T, "k:", lw=1)
    ax.set(xlabel=r"deflection $\theta$ [deg]", ylabel=r"wave angle $\beta$ [deg]")
    ax.legend()
    return ax


def plot_duct_flow(state: DuctFlowState, ax=None):
    """Plot the Fanno or Rayleigh property ratios against Mach number.

    Parameters
    ----------
    state : DuctFlowState
        From :func:`engineeringkit.aerodynamics.fanno_flow` or
        :func:`engineeringkit.aerodynamics.rayleigh_flow` with a Mach-number array.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    M = np.atleast_1d(state.mach)
    curves = [("$T/T^*$", state.temperature_ratio), ("$p/p^*$", state.pressure_ratio), ("$p_0/p_0^*$", state.stagnation_pressure_ratio)]
    curves.append(("$fL^*/D$", state.friction_parameter) if state.kind == "fanno" else ("$T_0/T_0^*$", state.stagnation_temperature_ratio))
    for label, values in curves:
        ax.plot(M, np.atleast_1d(values), label=label)
    ax.axvline(1, color="gray", lw=0.8)
    ax.set(xlabel="Mach number", ylabel="ratio to sonic state", title=f"{state.kind.capitalize()} flow")
    ax.legend()
    return ax


def plot_spanwise_loading(result: LiftingLineResult, ax=None):
    """Plot a lifting-line circulation distribution next to the elliptic load of equal lift.

    Parameters
    ----------
    result : LiftingLineResult
        From :func:`engineeringkit.aerodynamics.lifting_line`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    eta = 2 * result.y / result.span
    ax.plot(eta, result.circulation_per_speed, label=f"lifting line, e = {result.span_efficiency:.3f}")
    elliptic = 2 * result.span * result.coefficients[0] * np.sqrt(1 - eta**2)
    ax.plot(eta, elliptic, "--", label="elliptic, same lift")
    ax.set(xlabel="spanwise station $2y/b$", ylabel=r"$\Gamma / V_\infty$ [m]")
    ax.legend()
    return ax


def plot_pressure_coefficient(result: PanelResult, ax=None):
    """Plot :math:`-C_p` along the chord, upper surface solid and lower surface dashed.

    Parameters
    ----------
    result : PanelResult
        From :func:`engineeringkit.aerodynamics.panel_method`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    x = (result.x - result.nodes[:, 0].min()) / result.chord
    nose = int(np.argmin(result.nodes[:, 0]))
    ax.plot(x[nose:], -result.cp[nose:], label="upper surface")
    ax.plot(x[:nose], -result.cp[:nose], "--", label="lower surface")
    ax.set(xlabel="x / c", ylabel="$-C_p$", title=f"$c_l$ = {result.cl:.3f}")
    ax.legend()
    return ax
