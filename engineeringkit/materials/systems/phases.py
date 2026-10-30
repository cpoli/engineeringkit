r"""Phase diagrams and transformation kinetics: the lever rule and Avrami's law.

In a two-phase field of a binary phase diagram the tie line at the alloy's
temperature gives the compositions of the two phases, and a mass balance
gives their amounts, the lever rule. How fast a new phase forms is a
different question: Avrami (1939) showed that nuclei growing at random until
they impinge transform a fraction :math:`y = 1 - e^{-kt^n}`, the S-shaped
curve of recrystallization and of pearlite growing from austenite.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.materials.utils.fitting import power_law_fit

__all__ = ["lever_rule", "avrami_fraction", "avrami_time", "avrami_fit"]


def lever_rule(C0, C_alpha, C_beta):
    r"""Mass fractions of two coexisting phases from the tie line.

    .. math:: W_\alpha = \frac{C_\beta - C_0}{C_\beta - C_\alpha}, \qquad
              W_\beta = \frac{C_0 - C_\alpha}{C_\beta - C_\alpha}

    Parameters
    ----------
    C0 : float or array_like
        Overall alloy composition [wt fraction or wt%].
    C_alpha, C_beta : float or array_like
        Compositions of the two phases at the ends of the tie line [same
        unit], different; ``C0`` must lie between them.

    Returns
    -------
    tuple of float or ndarray
        ``(W_alpha, W_beta)``, mass fractions [-] summing to one.

    Notes
    -----
    Compositions in weight percent give mass fractions; atom percent gives
    mole fractions. Inputs are raw SI numbers with no dimensional metadata;
    see :doc:`/guides/materials`.

    Examples
    --------
    >>> W_L, W_alpha = lever_rule(35.0, 31.5, 42.5)
    >>> round(float(W_L), 2), round(float(W_alpha), 2)
    (0.68, 0.32)
    """
    _finite(C0=C0, C_alpha=C_alpha, C_beta=C_beta)
    C0, C_alpha, C_beta = np.asarray(C0), np.asarray(C_alpha), np.asarray(C_beta)
    if np.any(C_alpha == C_beta):
        raise ValueError("the tie line needs two different phase compositions")
    if np.any((C0 - C_alpha) * (C0 - C_beta) > 0):
        raise ValueError("C0 must lie on the tie line between C_alpha and C_beta")
    W_alpha = (C_beta - C0) / (C_beta - C_alpha)
    return W_alpha, 1 - W_alpha


def avrami_fraction(t, k, n):
    r"""Transformed fraction from the Avrami (Johnson-Mehl-Avrami-Kolmogorov) equation.

    .. math:: y = 1 - e^{-k t^n}

    Parameters
    ----------
    t : float or array_like
        Time since the transformation began [s], nonnegative.
    k : float or array_like
        Rate constant [1/s^n], positive.
    n : float or array_like
        Avrami exponent [-], positive (1 to 4: 4 for constant nucleation and
        three-dimensional growth, 3 for site-saturated growth).

    Returns
    -------
    float or ndarray
        Fraction transformed [-].

    Notes
    -----
    Isothermal transformation only. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    M. Avrami, "Kinetics of phase change. I. General theory," *J. Chem.
    Phys.* 7 (1939) 1103-1112.

    Examples
    --------
    >>> round(float(avrami_fraction(10.0, 1e-3, 3.0)), 4)
    0.6321
    """
    _nonnegative(t=t)
    _positive(k=k, n=n)
    return -np.expm1(-k * np.asarray(t, float) ** n)


def avrami_time(fraction, k, n):
    r"""Time to transform a given fraction, :math:`t = \left(-\ln(1 - y)/k\right)^{1/n}`.

    Parameters
    ----------
    fraction : float or array_like
        Fraction transformed :math:`y` [-], in (0, 1); 0.5 gives the half
        time whose inverse Callister calls the rate.
    k : float or array_like
        Rate constant [1/s^n], positive.
    n : float or array_like
        Avrami exponent [-], positive.

    Returns
    -------
    float or ndarray
        Time [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(avrami_time(0.5, 1e-3, 3.0)), 3)
    8.85
    """
    _interval("fraction", fraction, 0, 1, closed=False)
    _positive(k=k, n=n)
    return (-np.log1p(-np.asarray(fraction)) / k) ** (1 / np.asarray(n))


def avrami_fit(t, fraction):
    r"""Fit ``k`` and ``n`` to transformation data on Avrami axes.

    :math:`\ln(-\ln(1 - y)) = \ln k + n \ln t` is a straight line, so a
    power-law fit of :math:`-\ln(1 - y)` against :math:`t` gives both
    constants.

    Parameters
    ----------
    t : array_like
        Times [s], positive, at least two distinct.
    fraction : array_like
        Matching fractions transformed [-], in (0, 1).

    Returns
    -------
    tuple of float
        ``(k, n)``: the rate constant [1/s^n] and the exponent [-].

    Notes
    -----
    Points near 0 and 1 carry the most scatter; trim them before fitting.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _interval("fraction", fraction, 0, 1, closed=False)
    return power_law_fit(t, -np.log1p(-np.asarray(fraction, float)))
