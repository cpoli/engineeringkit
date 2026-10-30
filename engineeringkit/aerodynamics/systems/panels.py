r"""The Hess-Smith panel method for two-dimensional potential flow about a lifting body.

The surface is cut into straight panels, each carrying its own uniform source
strength :math:`q_j`, plus one vortex strength :math:`\gamma` shared by all
panels. Requiring zero normal velocity at every panel midpoint and equal
surface speeds on the two trailing-edge panels (the Kutta condition) gives
:math:`N + 1` linear equations for the :math:`N + 1` unknowns. Unlike thin-airfoil
theory the method keeps the real thickness, so it predicts the suction peak at
the leading edge and the extra lift a thick section carries.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit.aerodynamics.core.base import Airfoil, PanelResult

__all__ = ["panel_method"]


def _contour(geometry):
    """Panel end points, clockwise from the lower trailing edge, as an ``(N + 1, 2)`` array."""
    if isinstance(geometry, Airfoil):
        return np.vstack((geometry.lower[::-1], geometry.upper[1:]))
    nodes = np.asarray(geometry, float)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or len(nodes) < 4:
        raise ValueError("geometry must be an Airfoil or an (N + 1, 2) array of at least four points")
    _finite(geometry=nodes)
    return nodes


def _unit_panel_velocities(points, start, tangent, normal, length):
    """Velocity at each point from a unit source and a unit clockwise vortex on every panel.

    Returns two arrays of shape ``(points, panels, 2)``. A point on a panel
    (its own midpoint) takes the limit from the outward side.
    """
    d = points[:, None, :] - start[None, :, :]
    xp = np.einsum("ijk,jk->ij", d, tangent)
    yp = np.einsum("ijk,jk->ij", d, normal)
    own = np.isclose(yp, 0, atol=1e-12 * length.max()) & (xp > 0) & (xp < length)
    yp = np.where(own, 0.0, yp)
    u = np.log((xp**2 + yp**2) / ((xp - length) ** 2 + yp**2)) / (4 * np.pi)
    v = (np.arctan2(yp, xp - length) - np.arctan2(yp, xp)) / (2 * np.pi)
    v = np.where(own, 0.5, v)
    source = u[..., None] * tangent + v[..., None] * normal
    vortex = v[..., None] * tangent - u[..., None] * normal
    return source, vortex


def panel_method(geometry, alpha_rad) -> PanelResult:
    r"""Surface pressure, lift and moment of a 2-D body by the Hess-Smith source-vortex panel method.

    Parameters
    ----------
    geometry : Airfoil or array_like
        An :class:`~engineeringkit.aerodynamics.Airfoil` from
        :func:`~engineeringkit.aerodynamics.naca4_airfoil`, or the panel end
        points as an ``(N + 1, 2)`` array ordered clockwise from the lower
        surface at the trailing edge, round the leading edge, to the upper
        surface at the trailing edge. Repeat the first point last to close a
        smooth body such as a circle. Units are arbitrary but consistent.
    alpha_rad : float
        Angle of attack [rad].

    Returns
    -------
    PanelResult
        Surface speed and pressure coefficient at each panel midpoint, the
        singularity strengths, and the lift coefficient from the circulation
        and from the integrated pressure, with the quarter-chord moment.

    Raises
    ------
    ValueError
        If the points run counterclockwise or are not a list of 2-D points.

    Notes
    -----
    Steady, inviscid, incompressible flow: no boundary layer, separation or
    stall, and so no drag. The pressure-integrated lift converges to the
    circulation lift as the panels are refined. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/aerodynamics`.

    References
    ----------
    J. L. Hess and A. M. O. Smith, "Calculation of potential flow about
    arbitrary bodies," *Progress in Aerospace Sciences* 8 (1967) 1-138;
    J. D. Anderson, *Fundamentals of Aerodynamics*, 6th ed. (McGraw-Hill,
    2017), secs. 3.17 and 4.10.

    Examples
    --------
    >>> from engineeringkit.aerodynamics import naca4_airfoil
    >>> result = panel_method(naca4_airfoil("0012", 81), np.radians(5.0))
    >>> round(result.cl, 3)
    0.599
    """
    _finite(alpha_rad=alpha_rad)
    nodes = _contour(geometry)
    x, y = nodes[:, 0], nodes[:, 1]
    if np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]) + x[-1] * y[0] - x[0] * y[-1] >= 0:
        raise ValueError("geometry must run clockwise: lower surface first, from the trailing edge")
    start, delta = nodes[:-1], np.diff(nodes, axis=0)
    length = np.hypot(delta[:, 0], delta[:, 1])
    if np.any(length == 0):
        raise ValueError("geometry has repeated consecutive points")
    tangent = delta / length[:, None]
    normal = np.column_stack((-tangent[:, 1], tangent[:, 0]))
    mid = start + delta / 2

    source, vortex = _unit_panel_velocities(mid, start, tangent, normal, length)
    vortex = vortex.sum(axis=1)
    freestream = np.array([np.cos(alpha_rad), np.sin(alpha_rad)])
    n = len(length)
    matrix = np.zeros((n + 1, n + 1))
    rhs = np.zeros(n + 1)
    matrix[:n, :n] = np.einsum("ijk,ik->ij", source, normal)
    matrix[:n, n] = np.einsum("ik,ik->i", vortex, normal)
    rhs[:n] = -normal @ freestream
    tangential_source = np.einsum("ijk,ik->ij", source, tangent)
    tangential_vortex = np.einsum("ik,ik->i", vortex, tangent)
    matrix[n, :n] = tangential_source[0] + tangential_source[-1]
    matrix[n, n] = tangential_vortex[0] + tangential_vortex[-1]
    rhs[n] = -(tangent[0] + tangent[-1]) @ freestream
    solution = np.linalg.solve(matrix, rhs)
    q, gamma = solution[:n], solution[n]

    vt = tangent @ freestream + tangential_source @ q + gamma * tangential_vortex
    cp = 1 - vt**2
    chord = float(np.ptp(x))
    force = -(cp * length)[:, None] * normal
    lift_direction = np.array([-np.sin(alpha_rad), np.cos(alpha_rad)])
    leading_edge = nodes[int(np.argmin(x))]
    arm = mid - (leading_edge + [chord / 4, 0.0])
    moment = np.sum(arm[:, 0] * force[:, 1] - arm[:, 1] * force[:, 0])
    return PanelResult(
        nodes=nodes,
        x=mid[:, 0],
        y=mid[:, 1],
        tangential_velocity=vt,
        cp=cp,
        source_strengths=q,
        vortex_strength=float(gamma),
        alpha=float(alpha_rad),
        chord=chord,
        cl=float(2 * gamma * length.sum() / chord),
        cl_pressure=float(force.sum(axis=0) @ lift_direction / chord),
        cm_quarter_chord=float(-moment / chord**2),
    )
