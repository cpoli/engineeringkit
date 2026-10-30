r"""Static longitudinal stability: tail volume, downwash, neutral point and static margin.

An aircraft is statically stable in pitch when a nose-up disturbance creates
a nose-down moment, :math:`C_{M,\alpha} < 0`. With the wing-body lift acting
at its aerodynamic centre :math:`h_{ac,wb}` and the tail lift a distance
:math:`l_t` behind the centre of gravity, the moment slope is linear in the
centre-of-gravity position :math:`h = x_{cg}/\bar c` and vanishes at the
*neutral point* :math:`h_n`. The static margin :math:`h_n - h` measures how
far forward of it the centre of gravity sits.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = [
    "tail_volume_ratio",
    "downwash_gradient",
    "neutral_point",
    "pitching_moment_slope",
    "static_margin",
]


def tail_volume_ratio(S_t, l_t, S, c_bar):
    r"""Horizontal-tail volume ratio :math:`V_H = l_t S_t / (\bar c S)`.

    Parameters
    ----------
    S_t : float or array_like
        Horizontal-tail area [m^2], positive.
    l_t : float or array_like
        Distance from the centre of gravity to the tail aerodynamic centre [m], positive.
    S : float or array_like
        Wing area [m^2], positive.
    c_bar : float or array_like
        Mean aerodynamic chord [m], positive.

    Returns
    -------
    float or ndarray
        Tail volume ratio [-]; typically 0.3 to 1.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(tail_volume_ratio(3.0, 5.0, 16.0, 1.5)), 4)
    0.625
    """
    _positive(S_t=S_t, l_t=l_t, S=S, c_bar=c_bar)
    return np.asarray(l_t) * S_t / (c_bar * S)


def downwash_gradient(a_w, AR):
    r"""Downwash gradient at the tail, :math:`d\varepsilon/d\alpha = 2a_w/(\pi AR)`.

    Far behind an elliptically loaded wing the downwash angle is
    :math:`\varepsilon = 2C_L/(\pi AR)`, twice its value at the wing.

    Parameters
    ----------
    a_w : float or array_like
        Wing lift-curve slope [1/rad], positive.
    AR : float or array_like
        Wing aspect ratio [-], positive.

    Returns
    -------
    float or ndarray
        :math:`d\varepsilon/d\alpha` [-].

    Notes
    -----
    A far-wake estimate; a tail close behind the wing sees less. Inputs are
    raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(downwash_gradient(4.8, 8.0)), 4)
    0.382
    """
    _positive(a_w=a_w, AR=AR)
    return 2 * np.asarray(a_w) / (np.pi * AR)


def _check_tail(a_wb, a_t, V_H, deps_dalpha, eta):
    _positive(a_wb=a_wb, a_t=a_t, V_H=V_H)
    _interval("deps_dalpha", deps_dalpha, 0, 1, closed="left")
    _interval("eta", eta, 0, 1, closed="right")


def neutral_point(h_ac_wb, a_wb, a_t, V_H, deps_dalpha=0.0, eta=1.0):
    r"""Stick-fixed neutral point :math:`h_n = h_{ac,wb} + \eta V_H \frac{a_t}{a_{wb}}\left(1 - \frac{d\varepsilon}{d\alpha}\right)`.

    Parameters
    ----------
    h_ac_wb : float or array_like
        Wing-body aerodynamic centre as a fraction of the mean chord from
        its leading edge [-].
    a_wb : float or array_like
        Wing-body lift-curve slope [1/rad], positive.
    a_t : float or array_like
        Tail lift-curve slope [1/rad], positive.
    V_H : float or array_like
        Tail volume ratio [-], positive (:func:`tail_volume_ratio`).
    deps_dalpha : float or array_like, optional
        Downwash gradient at the tail [-] in ``[0, 1)`` (:func:`downwash_gradient`).
    eta : float or array_like, optional
        Tail efficiency, the tail-to-free-stream dynamic-pressure ratio [-], in ``(0, 1]``.

    Returns
    -------
    float or ndarray
        Neutral point as a fraction of the mean chord [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Introduction to Flight*, 8th ed. (McGraw-Hill, 2016),
    secs. 7.5-7.10.

    Examples
    --------
    >>> round(float(neutral_point(0.24, 4.8, 4.0, 0.6, 0.38)), 4)
    0.55
    """
    _finite(h_ac_wb=h_ac_wb)
    _check_tail(a_wb, a_t, V_H, deps_dalpha, eta)
    return np.asarray(h_ac_wb) + eta * np.asarray(V_H) * a_t / a_wb * (1 - np.asarray(deps_dalpha))


def pitching_moment_slope(h, h_ac_wb, a_wb, a_t, V_H, deps_dalpha=0.0, eta=1.0):
    r"""Moment-curve slope about the centre of gravity, :math:`C_{M,\alpha} = a_{wb}(h - h_n)`.

    Parameters
    ----------
    h : float or array_like
        Centre of gravity as a fraction of the mean chord [-].
    h_ac_wb, a_wb, a_t, V_H, deps_dalpha, eta
        As in :func:`neutral_point`.

    Returns
    -------
    float or ndarray
        :math:`dC_{M,cg}/d\alpha` [1/rad]; negative for a statically stable aircraft.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(pitching_moment_slope(0.3, 0.24, 4.8, 4.0, 0.6, 0.38)), 4)
    -1.2
    """
    _finite(h=h)
    return np.asarray(a_wb) * (np.asarray(h) - neutral_point(h_ac_wb, a_wb, a_t, V_H, deps_dalpha, eta))


def static_margin(h_n, h):
    r"""Static margin :math:`h_n - h`, the neutral point's distance behind the centre of gravity.

    Parameters
    ----------
    h_n : float or array_like
        Neutral point as a fraction of the mean chord [-].
    h : float or array_like
        Centre of gravity as a fraction of the mean chord [-].

    Returns
    -------
    float or ndarray
        Static margin [-]; positive for a statically stable aircraft
        (5-15% of the chord for most transport aircraft).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(static_margin(0.55, 0.3)), 4)
    0.25
    """
    _finite(h_n=h_n, h=h)
    return np.asarray(h_n) - h
