r"""First-order transients: the RC and RL time constants.

Kelvin's analysis of the first submarine telegraph cables showed that a
capacitance charged through a resistance responds on the time scale
:math:`RC`; an inductance in series with a resistance does the same on
:math:`L/R`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["rc_time_constant", "rl_time_constant", "rc_charging_voltage"]


def rc_time_constant(R, C):
    r"""RC time constant :math:`\tau = RC`.

    Parameters
    ----------
    R : float or array_like
        Resistance [ohm], positive.
    C : float or array_like
        Capacitance [F], positive.

    Returns
    -------
    float or ndarray
        Time constant [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    W. Thomson, "On the theory of the electric telegraph," *Proc. R. Soc.
    London* 7 (1855) 382-399.

    Examples
    --------
    >>> rc_time_constant(1e3, 1e-6)
    0.001
    """
    _finite(R=R, C=C)
    _positive(R=R, C=C)
    return R * C


def rl_time_constant(L, R):
    r"""RL time constant :math:`\tau = L/R`.

    Parameters
    ----------
    L : float or array_like
        Inductance [H], positive.
    R : float or array_like
        Resistance [ohm], positive.

    Returns
    -------
    float or ndarray
        Time constant [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.
    """
    _finite(L=L, R=R)
    _positive(L=L, R=R)
    return L / R


def rc_charging_voltage(v_source, R, C, t, v0=0.0):
    r"""Capacitor voltage :math:`v(t) = V_s + (v_0 - V_s)e^{-t/RC}` after a step.

    Parameters
    ----------
    v_source : float
        Constant source voltage applied at :math:`t = 0` [V].
    R : float
        Series resistance [ohm], positive.
    C : float
        Capacitance [F], positive.
    t : float or array_like
        Time since the step [s], nonnegative.
    v0 : float, optional
        Initial capacitor voltage [V].

    Returns
    -------
    ndarray
        Capacitor voltage [V], shaped like ``t``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> round(float(rc_charging_voltage(1.0, 1e3, 1e-6, 1e-3)), 4)
    0.6321
    """
    _finite(v_source=v_source, R=R, C=C, t=t, v0=v0)
    _nonnegative(t=t)
    t = np.asarray(t, float)
    tau = rc_time_constant(R, C)
    return v_source + (v0 - v_source) * np.exp(-t / tau)
