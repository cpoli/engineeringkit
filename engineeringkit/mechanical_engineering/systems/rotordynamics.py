r"""Rotordynamics: unbalance and the Jeffcott rotor.

An unbalanced disk on a flexible shaft whirls with radius
:math:`e r^2 / \sqrt{(1 - r^2)^2 + (2\zeta r)^2}`, where :math:`r` is the
spin speed over the critical speed :math:`\sqrt{k/m}`. Jeffcott showed that
above the critical speed the disk turns about its own centre of mass, so a
shaft can be run safely through and beyond its whirling speed.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.mechanical_engineering.core.base import JeffcottResult

__all__ = ["rotating_unbalance_force", "jeffcott_response"]


def rotating_unbalance_force(m_e, e, omega):
    r"""Centrifugal force :math:`m_e e \omega^2` of a rotating unbalance.

    Parameters
    ----------
    m_e : float or array_like
        Unbalanced mass [kg], nonnegative.
    e : float or array_like
        Its distance from the spin axis [m], nonnegative.
    omega : float or array_like
        Spin speed [rad/s].

    Returns
    -------
    float or ndarray
        Force amplitude [N], rotating with the shaft.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(rotating_unbalance_force(2.0, 0.1, 10.0)), 12)
    20.0
    """
    _finite(m_e=m_e, e=e, omega=omega)
    _nonnegative(m_e=m_e, e=e)
    return m_e * e * omega**2


def jeffcott_response(m, k, e, omega, c=0.0) -> JeffcottResult:
    r"""Steady synchronous whirl of a Jeffcott (Föppl-de Laval) rotor.

    A disk of mass :math:`m`, whose centre of mass sits :math:`e` from the
    shaft centre, spins at :math:`\omega` midway on a massless shaft of
    lateral stiffness :math:`k`. With :math:`r = \omega/\omega_c` and
    :math:`\zeta = c/(2\sqrt{km})`,

    .. math::

        \frac{|z|}{e} = \frac{r^2}{\sqrt{(1 - r^2)^2 + (2\zeta r)^2}}, \qquad
        \phi = \operatorname{atan2}(2\zeta r,\ 1 - r^2)

    Parameters
    ----------
    m : float
        Disk mass [kg], positive.
    k : float
        Shaft lateral stiffness at the disk [N/m], positive.
    e : float
        Mass eccentricity [m], positive.
    omega : float or array_like
        Spin speeds [rad/s], nonnegative.
    c : float, optional
        External viscous damping [N s/m], nonnegative.

    Returns
    -------
    JeffcottResult
        Whirl radius [m] and its lag behind the unbalance [rad] at each
        speed, and the critical speed [rad/s]. Without damping the radius
        is infinite at the critical speed.

    Notes
    -----
    ``m``, ``k``, ``e`` and ``c`` are scalar-only, since the result describes
    one rotor; ``omega`` broadcasts.

    Gyroscopic effects, internal damping and bearing flexibility are
    neglected. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    H. H. Jeffcott, "The lateral vibration of loaded shafts in the
    neighbourhood of a whirling speed -- the effect of want of balance,"
    *Phil. Mag.* (6) 37 (1919) 304-314.

    Examples
    --------
    >>> result = jeffcott_response(10.0, 1e6, 1e-4, [0.0, 1e4], c=100.0)
    >>> round(result.critical_speed, 6), round(float(result.amplitude[1] / 1e-4), 4)
    (316.227766, 1.001)
    """
    _positive(m=m, k=k, e=e)
    _nonnegative(omega=omega, c=c)
    w = np.asarray(omega, dtype=float)
    wc = float(np.sqrt(k / m))
    r = w / wc
    zeta = c / (2 * np.sqrt(k * m))
    with np.errstate(divide="ignore"):
        amplitude = e * r * r / np.hypot(1 - r * r, 2 * zeta * r)
    phase = np.arctan2(2 * zeta * r, 1 - r * r)
    return JeffcottResult(omega=w, amplitude=amplitude, phase=phase, critical_speed=wc, eccentricity=float(e))
