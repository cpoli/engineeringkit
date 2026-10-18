"""Plotting helpers for engineeringkit.structural_engineering: beam diagrams, trusses, frames and Mohr's circle."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering.systems.stress import mohr_circle

__all__ = ["plot_beam_diagrams", "plot_truss", "plot_frame", "plot_mohr_circle"]


def plot_beam_diagrams(result, axes=None):
    """Shear, moment and deflection diagrams of a :class:`BeamResult`, stacked.

    Parameters
    ----------
    result : BeamResult
        Output of :class:`~engineeringkit.structural_engineering.SimplySupportedBeam`.
    axes : sequence of three matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    ndarray of matplotlib.axes.Axes
    """
    if axes is None:
        _, axes = plt.subplots(3, 1, sharex=True)
    axes = np.asarray(axes)
    if axes.size != 3:
        raise ValueError("three axes are required")
    for ax, values in zip(axes, (result.shear, result.moment, result.deflection)):
        ax.plot(result.x, values)
        ax.axhline(0, color="gray", lw=0.5)
    axes[-1].set_xlabel("x [m]")
    return axes


def plot_truss(nodes, members, displacements=None, forces=None, scale=1.0, ax=None):
    """Draw a plane truss, optionally deformed and coloured by member force.

    Parameters
    ----------
    nodes : array_like, shape (n, 2)
        Joint coordinates [m].
    members : sequence of (int, int)
        Member connectivity.
    displacements : array_like, optional
        Joint displacements [m], two per joint. If given, the undeformed truss
        is dashed and the deformed one solid.
    forces : array_like, optional
        Member forces [N]; tension is drawn red and compression blue.
    scale : float, optional
        Displacement magnification [-].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    xy = np.asarray(nodes, float)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError("nodes must have shape (n, 2)")
    if ax is None:
        _, ax = plt.subplots()
    shown = xy
    if displacements is not None:
        shown = xy + scale * np.asarray(displacements, float).reshape(-1, 2)
    signs = np.zeros(len(members)) if forces is None else np.sign(np.asarray(forces, float))
    for (i, j), sign in zip(members, signs):
        if displacements is not None:
            ax.plot(xy[[i, j], 0], xy[[i, j], 1], "--", color="gray", alpha=0.5)
        color = "C3" if sign > 0 else "C0" if sign < 0 else "k"
        ax.plot(shown[[i, j], 0], shown[[i, j], 1], color=color)
    ax.set_aspect("equal")
    ax.set(xlabel="x [m]", ylabel="y [m]")
    return ax


def plot_frame(result, scale=1.0, ax=None, n=21):
    """Draw a plane frame and its deflected shape from a :class:`FrameResult`.

    Each member's deflected shape is the cubic Hermite interpolation of its
    end displacements and rotations, which is exact for members without
    span loads.

    Parameters
    ----------
    result : FrameResult
        Output of :meth:`PlanarFrame.solve`.
    scale : float, optional
        Displacement magnification [-].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    n : int, optional
        Points drawn along each member.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    xy, d = result.nodes, result.displacements.reshape(-1, 3)
    t = np.linspace(0.0, 1.0, n)
    h1, h2, h3, h4 = 1 - 3 * t**2 + 2 * t**3, t - 2 * t**2 + t**3, 3 * t**2 - 2 * t**3, -(t**2) + t**3
    for i, j in result.members:
        dx, dy = xy[j] - xy[i]
        L = np.hypot(dx, dy)
        c, s = dx / L, dy / L
        ui, vi = c * d[i, 0] + s * d[i, 1], -s * d[i, 0] + c * d[i, 1]
        uj, vj = c * d[j, 0] + s * d[j, 1], -s * d[j, 0] + c * d[j, 1]
        along = t * L + scale * ((1 - t) * ui + t * uj)
        across = scale * (h1 * vi + h2 * L * d[i, 2] + h3 * vj + h4 * L * d[j, 2])
        ax.plot(xy[[i, j], 0], xy[[i, j], 1], "--", color="gray", alpha=0.5)
        ax.plot(xy[i, 0] + c * along - s * across, xy[i, 1] + s * along + c * across, color="C0")
    ax.set_aspect("equal")
    ax.set(xlabel="x [m]", ylabel="y [m]")
    return ax


def plot_mohr_circle(sx, sy, txy, ax=None):
    """Draw Mohr's circle for a plane stress state.

    The x-face point :math:`(\\sigma_x, \\tau_{xy})` and the y-face point
    :math:`(\\sigma_y, -\\tau_{xy})` are marked, with shear plotted positive
    downward so that rotating the element counterclockwise by :math:`\\theta`
    moves the point counterclockwise by :math:`2\\theta`.

    Parameters
    ----------
    sx, sy, txy : float
        Plane stress components [Pa].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    center, radius, s1, s2 = (float(v) for v in mohr_circle(sx, sy, txy))
    if ax is None:
        _, ax = plt.subplots()
    phi = np.linspace(0, 2 * np.pi, 361)
    ax.plot(center + radius * np.cos(phi), radius * np.sin(phi), color="C0")
    ax.plot([sx, sy], [txy, -txy], "o-", color="C1", label="x and y faces")
    ax.plot([s1, s2], [0, 0], "s", color="C3", label="principal stresses")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_aspect("equal")
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    ax.set(xlabel="normal stress [Pa]", ylabel="shear stress [Pa] (positive down)")
    ax.legend()
    return ax
