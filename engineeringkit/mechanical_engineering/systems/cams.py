r"""Cam displacement laws: follower lift and its derivatives over one rise.

A cam turning through the rise angle :math:`\beta` lifts its follower by
:math:`L`. The law :math:`s(\theta)` fixes the follower's velocity,
acceleration and jerk, which are :math:`\omega`, :math:`\omega^2` and
:math:`\omega^3` times the derivatives with respect to cam angle. Uniform
velocity has infinite acceleration at the ends; parabolic (constant
acceleration) and simple harmonic motion have finite acceleration but
infinite jerk where they join a dwell; cycloidal motion and the 3-4-5 and
4-5-6-7 polynomials start and end with zero acceleration.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit.mechanical_engineering.core.base import CamMotion

__all__ = ["cam_follower_motion"]


def _uniform(u):
    zero = np.zeros_like(u)
    return u, np.ones_like(u), zero, zero


def _parabolic(u):
    first = u < 0.5
    s = np.where(first, 2 * u * u, 1 - 2 * (1 - u) ** 2)
    v = np.where(first, 4 * u, 4 * (1 - u))
    return s, v, np.where(first, 4.0, -4.0), np.zeros_like(u)


def _harmonic(u):
    p = np.pi * u
    return (1 - np.cos(p)) / 2, np.pi / 2 * np.sin(p), np.pi**2 / 2 * np.cos(p), -(np.pi**3) / 2 * np.sin(p)


def _cycloidal(u):
    p = 2 * np.pi * u
    return u - np.sin(p) / (2 * np.pi), 1 - np.cos(p), 2 * np.pi * np.sin(p), 4 * np.pi**2 * np.cos(p)


def _polynomial(coefficients):
    s = np.polynomial.Polynomial(coefficients)
    derivatives = (s, s.deriv(1), s.deriv(2), s.deriv(3))
    return lambda u: tuple(p(u) for p in derivatives)


_LAWS = {
    "uniform": _uniform,
    "parabolic": _parabolic,
    "harmonic": _harmonic,
    "cycloidal": _cycloidal,
    "polynomial_345": _polynomial([0, 0, 0, 10, -15, 6]),
    "polynomial_4567": _polynomial([0, 0, 0, 0, 35, -84, 70, -20]),
}


def cam_follower_motion(theta, lift, rise_angle, law="cycloidal") -> CamMotion:
    r"""Follower displacement, velocity, acceleration and jerk over a cam rise.

    With :math:`u = \theta/\beta`, the laws are

    ===================  ================================================
    law                  :math:`s/L`
    ===================  ================================================
    ``uniform``          :math:`u`
    ``parabolic``        :math:`2u^2` then :math:`1 - 2(1 - u)^2`
    ``harmonic``         :math:`(1 - \cos\pi u)/2`
    ``cycloidal``        :math:`u - \sin(2\pi u)/2\pi`
    ``polynomial_345``   :math:`10u^3 - 15u^4 + 6u^5`
    ``polynomial_4567``  :math:`35u^4 - 84u^5 + 70u^6 - 20u^7`
    ===================  ================================================

    Parameters
    ----------
    theta : float or array_like
        Cam angle from the start of the rise [rad], in [0, ``rise_angle``].
    lift : float
        Total rise :math:`L` [m], positive.
    rise_angle : float
        Cam angle :math:`\beta` over which the rise happens [rad], positive.
    law : str, optional
        One of the names above.

    Returns
    -------
    CamMotion
        :math:`s` [m] and :math:`ds/d\theta`, :math:`d^2s/d\theta^2`,
        :math:`d^3s/d\theta^3` [m/rad^k]. The return from :math:`L` to zero
        is :math:`L - s`.

    Notes
    -----
    The jumps in acceleration (parabolic) or velocity (uniform) where a law
    meets a dwell are not represented as impulses. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> motion = cam_follower_motion([0.0, np.pi / 4, np.pi / 2], 0.01, np.pi / 2)
    >>> np.round(motion.displacement, 6).tolist()
    [0.0, 0.005, 0.01]
    """
    _positive(lift=lift, rise_angle=rise_angle)
    _interval("theta", theta, 0, rise_angle)
    if law not in _LAWS:
        raise ValueError(f"law must be one of {sorted(_LAWS)}")
    t = np.asarray(theta, dtype=float)
    s, v, a, j = (np.asarray(x, dtype=float) for x in _LAWS[law](t / rise_angle))
    return CamMotion(
        theta=t,
        displacement=lift * s,
        velocity=lift * v / rise_angle,
        acceleration=lift * a / rise_angle**2,
        jerk=lift * j / rise_angle**3,
        law=law,
    )
