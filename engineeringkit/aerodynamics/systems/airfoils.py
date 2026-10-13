r"""NACA four-digit airfoil sections.

Jacobs, Ward and Pinkerton (NACA Report 460, 1933) defined a family of
sections by a thickness distribution laid perpendicular to a two-parabola mean
camber line. The four digits ``MPTT`` give the maximum camber ``M`` (% chord),
its position ``P`` (tenths of chord) and the thickness ``TT`` (% chord).
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit.aerodynamics.core.base import Airfoil

__all__ = ["naca4_airfoil"]


def naca4_airfoil(designation: str, n: int = 101, closed_trailing_edge: bool = False) -> Airfoil:
    r"""Coordinates of a NACA four-digit section on a unit chord.

    The half-thickness is

    .. math::

        y_t = 5t\,(0.2969\sqrt{x} - 0.1260x - 0.3516x^2 + 0.2843x^3 - 0.1015x^4),

    and the surfaces are offset from the camber line along its normal.
    Stations are cosine-spaced to resolve the leading edge.

    Parameters
    ----------
    designation : str
        Four digits, e.g. ``"0012"`` or ``"2412"``.
    n : int, optional
        Number of chordwise stations, at least 3.
    closed_trailing_edge : bool, optional
        Use ``-0.1036`` instead of ``-0.1015`` for the quartic coefficient so
        the trailing edge has zero thickness.

    Returns
    -------
    Airfoil

    Notes
    -----
    Coordinates are fractions of chord; see :doc:`/guides/aerodynamics`.

    References
    ----------
    E. N. Jacobs, K. E. Ward and R. M. Pinkerton, "The characteristics of 78
    related airfoil sections from tests in the variable-density wind tunnel,"
    NACA Report 460 (1933).

    Examples
    --------
    >>> foil = naca4_airfoil("2412")
    >>> round(foil.summary()["max_thickness"], 3)
    0.12
    """
    if not isinstance(designation, str) or len(designation) != 4 or not designation.isdigit():
        raise ValueError("designation must be four digits, e.g. '2412'")
    _integer("n", n, 3)
    m, p, t = int(designation[0]) / 100, int(designation[1]) / 10, int(designation[2:]) / 100
    if m > 0:
        _interval("camber position", p, 0, 1, closed=False)
    if t == 0:
        raise ValueError("thickness digits must be nonzero")
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))
    a4 = -0.1036 if closed_trailing_edge else -0.1015
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 + a4 * x**4)
    yc = np.zeros_like(x)
    dyc = np.zeros_like(x)
    if m > 0:
        front = x < p
        yc[front] = m / p**2 * (2 * p * x[front] - x[front] ** 2)
        dyc[front] = 2 * m / p**2 * (p - x[front])
        back = ~front
        yc[back] = m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x[back] - x[back] ** 2)
        dyc[back] = 2 * m / (1 - p) ** 2 * (p - x[back])
    theta = np.arctan(dyc)
    upper = np.column_stack((x - yt * np.sin(theta), yc + yt * np.cos(theta)))
    lower = np.column_stack((x + yt * np.sin(theta), yc - yt * np.cos(theta)))
    return Airfoil(f"NACA {designation}", x, yc, yt, upper, lower)
