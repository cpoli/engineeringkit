r"""Eigenvalues of the plane wall with a convective surface.

Separating variables in the transient plane wall of half-thickness :math:`L`
cooled by convection gives modes :math:`\cos(\zeta_n x/L)` whose
eigenvalues solve :math:`\zeta \tan\zeta = Bi`. Exactly one root lies in
each interval :math:`[(n-1)\pi, (n-1)\pi + \pi/2)`, so Brent's method on
:math:`\zeta\sin\zeta - Bi\cos\zeta` (continuous, unlike the tangent form)
brackets every root.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive

__all__ = ["plane_wall_eigenvalues"]


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
