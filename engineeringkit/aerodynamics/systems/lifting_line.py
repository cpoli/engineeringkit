r"""Prandtl's lifting-line equation for an arbitrary planform, solved by Glauert's Fourier series.

With the spanwise station :math:`y = -\tfrac b2\cos\theta` and the
circulation :math:`\Gamma(\theta) = 2bV_\infty\sum_n A_n\sin n\theta`, the
monoplane equation

.. math::

    \alpha(\theta) - \alpha_{L=0}(\theta) = \frac{4b}{a_0 c(\theta)}\sum_n A_n \sin n\theta
    + \sum_n n A_n \frac{\sin n\theta}{\sin\theta}

becomes a linear system when it is enforced at as many stations as there
are coefficients. The wing coefficients then follow from :math:`A_1` alone
(lift) and from all of them (induced drag).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive
from engineeringkit.aerodynamics.core.base import LiftingLineResult

__all__ = ["lifting_line"]


def _spanwise(name, value, y):
    """Evaluate a constant or a callable of the spanwise station ``y`` [m]."""
    values = np.broadcast_to(np.asarray(value(y) if callable(value) else value, float), y.shape)
    _finite(**{name: values})
    return values


def lifting_line(span, chord, alpha, alpha_l0=0.0, a0=2 * np.pi, n_terms=40) -> LiftingLineResult:
    r"""Spanwise lift distribution, lift and induced drag of an unswept wing.

    Parameters
    ----------
    span : float
        Wing span :math:`b` [m], positive.
    chord : float or callable
        Chord [m], positive: a constant, or a function of the spanwise
        station ``y`` [m] in :math:`[-b/2, b/2]` that accepts arrays.
    alpha : float or callable
        Geometric angle of attack [rad], a constant or a function of ``y``
        (washout makes it fall towards the tips).
    alpha_l0 : float or callable, optional
        Section zero-lift angle [rad], a constant or a function of ``y``.
    a0 : float, optional
        Section lift-curve slope [1/rad], positive; :math:`2\pi` from
        thin-airfoil theory.
    n_terms : int, optional
        Number of Fourier terms and collocation stations, at least one.

    Returns
    -------
    LiftingLineResult
        Fourier coefficients, the section lift and downwash at each
        collocation station, and the wing :math:`C_L`, :math:`C_{D,i}` and
        span efficiency.

    Notes
    -----
    Attached, incompressible flow over a straight wing of moderate to large
    aspect ratio; sweep, stall and the root and tip corrections of
    lifting-surface theory are not modelled. All harmonics are kept, so
    asymmetric twist (aileron deflection) is allowed. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/aerodynamics`.

    References
    ----------
    H. Glauert, *The Elements of Aerofoil and Airscrew Theory* (Cambridge
    University Press, 1926), ch. 11; J. D. Anderson, *Fundamentals of
    Aerodynamics*, 6th ed. (McGraw-Hill, 2017), sec. 5.3.

    Examples
    --------
    >>> wing = lifting_line(10.0, lambda y: 1.25 * np.sqrt(1 - (2 * y / 10.0) ** 2), 0.1)
    >>> round(wing.span_efficiency, 6), round(wing.CL, 5)
    (1.0, 0.5252)
    """
    _positive(span=span, a0=a0)
    _integer("n_terms", n_terms, minimum=1)
    n = np.arange(1, n_terms + 1)
    theta = np.arange(1, n_terms + 1) * np.pi / (n_terms + 1)
    y = -0.5 * span * np.cos(theta)
    c = _spanwise("chord", chord, y)
    _positive(chord=c)
    rhs = (_spanwise("alpha", alpha, y) - _spanwise("alpha_l0", alpha_l0, y)) * np.sin(theta)
    sines = np.sin(np.outer(theta, n))
    matrix = sines * (4 * span * np.sin(theta) / (a0 * c))[:, None] + sines * n
    A = np.linalg.solve(matrix, rhs)

    def chord_at(t):
        return float(np.asarray(chord(-0.5 * span * np.cos(t)) if callable(chord) else chord))

    area = quad(lambda t: chord_at(t) * np.sin(t), 0, np.pi)[0] * span / 2
    AR = span**2 / area
    weighted = float(np.sum(n * A**2))
    gamma_series = sines @ A
    return LiftingLineResult(
        span=float(span),
        area=float(area),
        aspect_ratio=float(AR),
        coefficients=A,
        y=y,
        chord=c,
        cl=4 * span * gamma_series / c,
        induced_angle=(sines * n) @ A / np.sin(theta),
        circulation_per_speed=2 * span * gamma_series,
        CL=float(np.pi * AR * A[0]),
        CDi=float(np.pi * AR * weighted),
        span_efficiency=float(A[0] ** 2 / weighted) if weighted else float("nan"),
    )
