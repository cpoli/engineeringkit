r"""Eigenvalues of the plane wall, cylinder and sphere with a convective surface.

Separating variables in the transient plane wall of half-thickness :math:`L`
cooled by convection gives modes :math:`\cos(\zeta_n x/L)` whose
eigenvalues solve :math:`\zeta \tan\zeta = Bi`. Exactly one root lies in
each interval :math:`[(n-1)\pi, (n-1)\pi + \pi/2)`, so Brent's method on
:math:`\zeta\sin\zeta - Bi\cos\zeta` (continuous, unlike the tangent form)
brackets every root. The cylinder's :math:`\zeta J_1(\zeta)/J_0(\zeta) = Bi`
has one root between each zero of :math:`J_1` and the next zero of
:math:`J_0`, and the sphere's :math:`1 - \zeta\cot\zeta = Bi` one root in
each :math:`((n-1)\pi, n\pi)`.
"""

from __future__ import annotations

import numpy as np
from scipy import special
from scipy.optimize import brentq

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive

__all__ = ["plane_wall_eigenvalues", "cylinder_eigenvalues", "sphere_eigenvalues"]


def plane_wall_eigenvalues(Bi, terms):
    r"""First ``terms`` positive roots of :math:`\zeta \tan\zeta = Bi`.

    Parameters
    ----------
    Bi : float
        Biot number :math:`hL/k` [-] on the wall half-thickness, positive.
    terms : int
        Number of eigenvalues, at least one.

    Returns
    -------
    ndarray
        Eigenvalues :math:`\zeta_1 < \zeta_2 < \cdots` [-]; :math:`\zeta_n`
        lies in :math:`((n-1)\pi, (n-1)\pi + \pi/2)`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(plane_wall_eigenvalues(1.0, 1)[0]), 4)
    0.8603
    """
    _positive(Bi=Bi)
    _integer("terms", terms, 1)
    Bi = float(Bi)

    def f(z):
        return z * np.sin(z) - Bi * np.cos(z)

    roots = [brentq(f, n * np.pi, n * np.pi + np.pi / 2, xtol=1e-14) for n in range(terms)]
    return np.array(roots)


def cylinder_eigenvalues(Bi, terms):
    r"""First ``terms`` positive roots of :math:`\zeta J_1(\zeta) = Bi\,J_0(\zeta)`.

    Parameters
    ----------
    Bi : float
        Biot number :math:`hr_o/k` [-] on the outer radius, positive.
    terms : int
        Number of eigenvalues, at least one.

    Returns
    -------
    ndarray
        Eigenvalues :math:`\zeta_1 < \zeta_2 < \cdots` [-]; :math:`\zeta_n`
        lies between the :math:`(n-1)`-th zero of :math:`J_1` (0 for
        :math:`n = 1`) and the :math:`n`-th zero of :math:`J_0`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(cylinder_eigenvalues(1.0, 1)[0]), 4)
    1.2558
    """
    _positive(Bi=Bi)
    _integer("terms", terms, 1)
    Bi = float(Bi)

    def f(z):
        return z * special.j1(z) - Bi * special.j0(z)

    lower = np.concatenate([[0.0], special.jn_zeros(1, terms - 1)]) if terms > 1 else np.zeros(1)
    upper = special.jn_zeros(0, terms)
    return np.array([brentq(f, a, b, xtol=1e-14) for a, b in zip(lower, upper)])


def sphere_eigenvalues(Bi, terms):
    r"""First ``terms`` positive roots of :math:`1 - \zeta\cot\zeta = Bi`.

    Parameters
    ----------
    Bi : float
        Biot number :math:`hr_o/k` [-] on the outer radius, positive.
    terms : int
        Number of eigenvalues, at least one.

    Returns
    -------
    ndarray
        Eigenvalues :math:`\zeta_1 < \zeta_2 < \cdots` [-]; :math:`\zeta_n`
        lies in :math:`((n-1)\pi, n\pi)`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Solved as
    :math:`\zeta\cos\zeta + (Bi - 1)\sin\zeta = 0`, which is positive
    just above zero (where it behaves as :math:`Bi\,\zeta`).

    Examples
    --------
    >>> round(float(sphere_eigenvalues(1.0, 1)[0]), 4)
    1.5708
    """
    _positive(Bi=Bi)
    _integer("terms", terms, 1)
    Bi = float(Bi)

    def f(z):
        return z * np.cos(z) + (Bi - 1) * np.sin(z)

    # f(0) = 0, so the first bracket starts just above zero but below the
    # first root, which is about sqrt(3 Bi) for small Bi.
    lower = np.arange(terms) * np.pi
    lower[0] = 1e-3 * min(1.0, np.sqrt(Bi))
    return np.array([brentq(f, a, (n + 1) * np.pi, xtol=1e-14) for n, a in enumerate(lower)])
