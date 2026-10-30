r"""Helical compression springs: rate, index and the Wahl stress factor.

A close-coiled helical spring is a wire in torsion: the axial load
:math:`F` twists it with torque :math:`FD/2`, so the rate is
:math:`k = Gd^4/8D^3N` and the nominal shear stress :math:`8FD/\pi d^3`.
The direct shear of the load and the curvature of the coil raise the stress
on the inside of the coil; Wahl's factor :math:`(4C - 1)/(4C - 4) + 0.615/C`
corrects for both, with :math:`C = D/d` the spring index.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = ["spring_index", "wahl_factor", "bergstrasser_factor", "spring_rate", "spring_shear_stress"]


def _index(C):
    _interval("C", C, 1, np.inf, closed=False)
    return np.asarray(C, dtype=float)


def spring_index(D, d):
    r"""Spring index :math:`C = D/d`, mean coil diameter over wire diameter.

    Parameters
    ----------
    D : float or array_like
        Mean coil diameter [m], positive.
    d : float or array_like
        Wire diameter [m], positive and less than ``D``.

    Returns
    -------
    float or ndarray
        :math:`C` [-]; practical springs have :math:`4 \le C \le 12`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> float(spring_index(0.03, 0.005))
    6.0
    """
    _positive(D=D, d=d)
    C = np.asarray(D, dtype=float) / d
    if np.any(C <= 1):
        raise ValueError("the wire diameter d must be smaller than the coil diameter D")
    return C


def wahl_factor(C):
    r"""Wahl's stress-correction factor :math:`K_W = \frac{4C - 1}{4C - 4} + \frac{0.615}{C}`.

    Parameters
    ----------
    C : float or array_like
        Spring index [-], greater than one.

    Returns
    -------
    float or ndarray
        :math:`K_W` [-], the inner-fibre shear over :math:`8FD/\pi d^3`,
        including direct shear and coil curvature.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    A. M. Wahl, "Stresses in heavy closely coiled helical springs," *Trans.
    ASME* 51 (1929) APM-51-17; A. M. Wahl, *Mechanical Springs* (Penton,
    Cleveland, 1944).

    Examples
    --------
    >>> round(float(wahl_factor(6.0)), 4)
    1.2525
    """
    c = _index(C)
    return (4 * c - 1) / (4 * c - 4) + 0.615 / c


def bergstrasser_factor(C):
    r"""Bergsträsser's factor :math:`K_B = (4C + 2)/(4C - 3)`, a simpler fit to the same stress.

    Parameters
    ----------
    C : float or array_like
        Spring index [-], greater than one.

    Returns
    -------
    float or ndarray
        :math:`K_B` [-], 0.7% to 1.5% below :func:`wahl_factor` for
        :math:`3 \le C \le 12`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(bergstrasser_factor(6.0)), 4)
    1.2381
    """
    c = _index(C)
    return (4 * c + 2) / (4 * c - 3)


def spring_rate(d, D, active_coils, G):
    r"""Rate of a close-coiled helical spring, :math:`k = Gd^4/(8D^3N_a)`.

    Parameters
    ----------
    d : float or array_like
        Wire diameter [m], positive.
    D : float or array_like
        Mean coil diameter [m], positive.
    active_coils : float or array_like
        Number of active coils :math:`N_a` [-], positive.
    G : float or array_like
        Shear modulus of the wire [Pa], positive.

    Returns
    -------
    float or ndarray
        Rate [N/m].

    Notes
    -----
    Only torsion of the wire is counted; the direct-shear deflection, a
    fraction of about :math:`1/2C^2`, is neglected. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(spring_rate(0.005, 0.03, 10, 79.3e9)), 3)
    22945.602
    """
    _positive(d=d, D=D, active_coils=active_coils, G=G)
    return G * np.asarray(d, dtype=float) ** 4 / (8 * np.asarray(D, dtype=float) ** 3 * active_coils)


def spring_shear_stress(F, D, d, correction="wahl"):
    r"""Maximum shear stress in the wire of a helical spring, :math:`\tau = K\,8FD/\pi d^3`.

    Parameters
    ----------
    F : float or array_like
        Axial load [N].
    D : float or array_like
        Mean coil diameter [m], positive.
    d : float or array_like
        Wire diameter [m], positive and less than ``D``.
    correction : {"wahl", "bergstrasser", "direct", "none"}, optional
        Factor :math:`K`: Wahl's or Bergsträsser's (direct shear and
        curvature), direct shear only, :math:`K_s = 1 + 0.5/C`, or none.

    Returns
    -------
    float or ndarray
        Shear stress [Pa], with the sign of ``F``.

    Notes
    -----
    Use the curvature-corrected factors for fatigue; for static loads of
    ductile wire, where local yielding relieves the curvature peak, the
    direct-shear factor is usual. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(spring_shear_stress(100.0, 0.03, 0.005, correction="none")) / 1e6, 4)
    61.1155
    """
    _finite(F=F)
    C = spring_index(D, d)
    factors = {
        "wahl": wahl_factor,
        "bergstrasser": bergstrasser_factor,
        "direct": lambda c: 1 + 0.5 / c,
        "none": lambda c: np.ones_like(c),
    }
    if correction not in factors:
        raise ValueError(f"correction must be one of {sorted(factors)}")
    return factors[correction](C) * 8 * np.asarray(F, dtype=float) * D / (np.pi * np.asarray(d, dtype=float) ** 3)
