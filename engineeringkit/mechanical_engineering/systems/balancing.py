r"""Balancing of rotating masses in one and two correction planes.

A mass :math:`m_i` at radius :math:`r_i` and angle :math:`\theta_i` on a
shaft spinning at :math:`\omega` pulls on the bearings with the rotating
force :math:`m_i r_i\omega^2`. Static balance makes the vector sum
:math:`\sum m_i r_i e^{i\theta_i}` vanish, so the shaft rests in any
position. Masses spread along the shaft also produce a rocking couple;
dynamic balance cancels both force and couple, which takes corrections in
two planes. Dalby's reference-plane method takes moments about one
correction plane to find the other.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import vector as _vector
from engineeringkit.mechanical_engineering.core.base import BalanceResult

__all__ = ["unbalance", "balance_rotating_masses"]


def _masses(masses, radii, angles):
    m, r, a = _vector("masses", masses), _vector("radii", radii), _vector("angles", angles)
    if not m.size == r.size == a.size:
        raise ValueError("masses, radii and angles must have the same length")
    _nonnegative(masses=m, radii=r)
    return m * r * np.exp(1j * a)


def unbalance(masses, radii, angles, positions=None):
    r"""Resultant unbalance and, optionally, its moment about the axial origin.

    Parameters
    ----------
    masses : array_like
        Rotating masses [kg], nonnegative.
    radii : array_like
        Their distances from the spin axis [m], nonnegative.
    angles : array_like
        Their angular positions [rad].
    positions : array_like, optional
        Their axial positions along the shaft [m].

    Returns
    -------
    float or tuple of float
        Magnitude of :math:`\sum m_i r_i e^{i\theta_i}` [kg m]; with
        ``positions``, also the magnitude of the couple unbalance
        :math:`\sum m_i r_i z_i e^{i\theta_i}` [kg m^2]. Multiply by
        :math:`\omega^2` for the shaking force [N] and moment [N m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(unbalance([1.0, 1.0], [0.1, 0.1], [0.0, np.pi]), 12)
    0.0
    """
    u = _masses(masses, radii, angles)
    if positions is None:
        return float(abs(u.sum()))
    z = _vector("positions", positions)
    if z.size != u.size:
        raise ValueError("positions must have one entry per mass")
    return float(abs(u.sum())), float(abs((u * z).sum()))


def balance_rotating_masses(masses, radii, angles, positions=None, planes=None) -> BalanceResult:
    r"""Correction masses for static (one-plane) or dynamic (two-plane) balance.

    Without ``positions`` and ``planes`` one correction cancels the force,
    :math:`m_c r_c e^{i\theta_c} = -\sum m_i r_i e^{i\theta_i}`. With them,
    corrections :math:`U_L` and :math:`U_R` at :math:`z_L` and :math:`z_R`
    cancel force and couple: moments about plane :math:`L` give

    .. math::

        U_R = -\frac{\sum m_i r_i (z_i - z_L) e^{i\theta_i}}{z_R - z_L},
        \qquad U_L = -\sum m_i r_i e^{i\theta_i} - U_R

    Parameters
    ----------
    masses : array_like
        Rotating masses [kg], nonnegative.
    radii : array_like
        Their radii [m], nonnegative.
    angles : array_like
        Their angular positions [rad].
    positions : array_like, optional
        Their axial positions [m]; give with ``planes`` for dynamic balance.
    planes : array_like, optional
        Axial positions :math:`(z_L, z_R)` of the two correction planes [m],
        distinct.

    Returns
    -------
    BalanceResult
        Correction mass-radius product [kg m] and angle [rad] in each plane.

    Notes
    -----
    The rotor is rigid: the corrections hold at every speed below the first
    bending critical speed. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/mechanical_engineering`.

    References
    ----------
    W. E. Dalby, *The Balancing of Engines* (Edward Arnold, London, 1902).

    Examples
    --------
    >>> result = balance_rotating_masses([2.0], [0.1], [np.pi / 2])
    >>> round(float(result.mass_radius[0]), 12), round(float(result.angle[0]), 12)
    (0.2, -1.570796326795)
    """
    u = _masses(masses, radii, angles)
    if (positions is None) != (planes is None):
        raise ValueError("give positions and planes together for two-plane balance")
    if positions is None:
        corrections = np.array([-u.sum()])
        plane_positions = np.array([])
    else:
        z = _vector("positions", positions)
        plane_positions = _vector("planes", planes)
        if z.size != u.size or plane_positions.size != 2:
            raise ValueError("positions needs one entry per mass and planes exactly two")
        zl, zr = plane_positions
        _finite(planes=plane_positions)
        if zl == zr:
            raise ValueError("the two correction planes must be distinct")
        right = -(u * (z - zl)).sum() / (zr - zl)
        corrections = np.array([-u.sum() - right, right])
    return BalanceResult(mass_radius=np.abs(corrections), angle=np.angle(corrections), plane_positions=plane_positions)
