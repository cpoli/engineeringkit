r"""Pressure vessels: thin-walled cylinders and spheres, and Lamé's thick-walled cylinder.

Cutting a thin cylinder of radius :math:`r` and wall :math:`t` along its
length and across it gives the hoop stress :math:`pr/t` and the
longitudinal stress :math:`pr/(2t)`; a sphere carries :math:`pr/(2t)` in
every direction. In a thick wall the hoop stress is no longer uniform:
Lamé's solution of the axisymmetric elasticity problem gives
:math:`\sigma_r = A - B/r^2` and :math:`\sigma_\theta = A + B/r^2`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["thin_wall_vessel_stresses", "thick_wall_cylinder_stresses"]

_SHAPES = ("cylinder", "sphere")


def thin_wall_vessel_stresses(p, r, t, shape="cylinder"):
    r"""Membrane stresses in a thin-walled pressure vessel.

    Parameters
    ----------
    p : float or array_like
        Internal gauge pressure [Pa].
    r : float or array_like
        Inner radius [m], positive.
    t : float or array_like
        Wall thickness [m], positive.
    shape : {"cylinder", "sphere"}, optional
        Vessel shape.

    Returns
    -------
    tuple of float or ndarray
        ``(hoop, longitudinal)`` stresses [Pa]: :math:`(pr/t, pr/2t)` for a
        cylinder with closed ends and :math:`(pr/2t, pr/2t)` for a sphere.
        The radial stress, between :math:`-p` and zero, is neglected.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. The usual limit is :math:`r/t \ge 10`, where the
    error against :func:`thick_wall_cylinder_stresses` is about 5%.
    """
    _finite(p=p)
    _positive(r=r, t=t)
    if shape not in _SHAPES:
        raise ValueError(f"shape must be one of {_SHAPES}")
    membrane = np.asarray(p) * r / (2 * np.asarray(t))
    return (membrane if shape == "sphere" else 2 * membrane), membrane


def thick_wall_cylinder_stresses(r, inner_radius, outer_radius, inner_pressure, outer_pressure=0.0):
    r"""Lamé's radial, hoop and axial stresses in a thick-walled cylinder.

    .. math::

        \sigma_r = A - \frac{B}{r^2},\qquad \sigma_\theta = A + \frac{B}{r^2},\qquad
        A = \frac{p_ir_i^2 - p_or_o^2}{r_o^2 - r_i^2},\quad
        B = \frac{(p_i - p_o)r_i^2r_o^2}{r_o^2 - r_i^2},

    and, with closed ends, :math:`\sigma_z = A`.

    Parameters
    ----------
    r : float or array_like
        Radius [m], in ``[inner_radius, outer_radius]``.
    inner_radius : float
        Inner radius :math:`r_i` [m], positive.
    outer_radius : float
        Outer radius :math:`r_o` [m], greater than ``inner_radius``.
    inner_pressure : float
        Internal pressure :math:`p_i` [Pa], nonnegative.
    outer_pressure : float, optional
        External pressure :math:`p_o` [Pa], nonnegative.

    Returns
    -------
    tuple of float or ndarray
        ``(radial, hoop, axial)`` stresses [Pa], tension positive. The hoop
        stress is largest at the bore.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Linear elastic, isotropic, away from the ends.

    References
    ----------
    G. Lamé and É. Clapeyron, "Mémoire sur l'équilibre intérieur des corps
    solides homogènes," *Mémoires présentés par divers savans à l'Académie
    royale des sciences de l'Institut de France* 4 (1833) 465-562.
    """
    _positive(inner_radius=inner_radius)
    _positive(wall=outer_radius - inner_radius)
    _nonnegative(inner_pressure=inner_pressure, outer_pressure=outer_pressure)
    _interval("r", r, inner_radius, outer_radius)
    ri2, ro2 = inner_radius**2, outer_radius**2
    A = (inner_pressure * ri2 - outer_pressure * ro2) / (ro2 - ri2)
    B = (inner_pressure - outer_pressure) * ri2 * ro2 / (ro2 - ri2)
    r2 = np.asarray(r, float) ** 2
    return A - B / r2, A + B / r2, A + 0 * r2
