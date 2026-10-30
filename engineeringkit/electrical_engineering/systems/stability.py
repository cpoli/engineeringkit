r"""Transient stability: the swing equation and the equal-area criterion.

A synchronous generator feeding an infinite bus through a reactance
delivers :math:`P_e = P_{max}\sin\delta`. A fault lowers :math:`P_{max}`,
the rotor accelerates under the unbalanced mechanical power, and once the
fault is cleared it must give back the kinetic energy it gained before its
angle passes the unstable equilibrium. Integrating the swing equation gives
the swing curve; comparing the accelerating and decelerating areas under
the power-angle curves gives the critical clearing angle without it.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.electrical_engineering.core.base import SwingResult
from engineeringkit.integrators import integrate_ode

__all__ = ["equal_area_clearing_angle", "critical_clearing_time", "swing_response"]


def _check_powers(p_mech, p_max_pre, p_max_fault, p_max_post):
    _finite(p_mech=p_mech, p_max_pre=p_max_pre, p_max_fault=p_max_fault, p_max_post=p_max_post)
    _positive(p_mech=p_mech)
    _nonnegative(p_max_fault=p_max_fault)
    if np.any(np.asarray(p_max_pre) <= p_mech) or np.any(np.asarray(p_max_post) <= p_mech):
        raise ValueError("p_max_pre and p_max_post must exceed p_mech, or there is no stable operating point")
    if np.any(np.asarray(p_max_fault) >= p_max_post):
        raise ValueError("p_max_fault must be below p_max_post")


def equal_area_clearing_angle(p_mech, p_max_pre, p_max_fault, p_max_post):
    r"""Initial and critical clearing rotor angles by the equal-area criterion.

    The rotor starts at :math:`\delta_0 = \arcsin(P_m/P_1)` on the prefault
    curve, accelerates along the fault curve :math:`P_2\sin\delta` and, once
    cleared, decelerates along :math:`P_3\sin\delta` until at most
    :math:`\delta_{max} = \pi - \arcsin(P_m/P_3)`. Equating the accelerating
    and decelerating areas,

    .. math::

        \cos\delta_{cr} = \frac{P_m(\delta_{max} - \delta_0)
        + P_3\cos\delta_{max} - P_2\cos\delta_0}{P_3 - P_2}.

    Parameters
    ----------
    p_mech : float or array_like
        Mechanical input power :math:`P_m` [pu], positive.
    p_max_pre, p_max_fault, p_max_post : float or array_like
        Amplitudes :math:`P_1, P_2, P_3` of the prefault, fault-on and
        post-fault power-angle curves [pu]; :math:`P_1, P_3 > P_m` and
        :math:`0 \le P_2 < P_3`.

    Returns
    -------
    tuple of float or ndarray
        ``(delta0, delta_critical)`` [rad]. If the fault-on curve keeps the
        rotor from passing :math:`\delta_{max}` whenever it is cleared,
        ``delta_critical`` is :math:`\delta_{max}`.

    Raises
    ------
    ValueError
        For inconsistent powers, or if the post-fault system cannot hold the
        rotor even when the fault is cleared at once.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. Classical model: constant
    voltage behind transient reactance, constant mechanical power, no
    damping, one machine against an infinite bus.

    References
    ----------
    R. D. Evans and C. F. Wagner, "Studies of transmission stability,"
    *Trans. AIEE* 45 (1926); E. W. Kimbark, *Power System Stability*,
    vol. 1 (New York: Wiley, 1948), ch. 4.

    Examples
    --------
    >>> d0, dcr = equal_area_clearing_angle(1.0, 2.0, 0.0, 2.0)
    >>> round(float(np.degrees(d0)), 2), round(float(np.degrees(dcr)), 2)
    (30.0, 79.56)
    """
    _check_powers(p_mech, p_max_pre, p_max_fault, p_max_post)
    p_mech, p_max_pre, p_max_fault, p_max_post = (np.asarray(x, dtype=float) for x in (p_mech, p_max_pre, p_max_fault, p_max_post))
    delta0 = np.arcsin(p_mech / p_max_pre)
    delta_max = np.pi - np.arcsin(p_mech / p_max_post)
    cos_cr = (p_mech * (delta_max - delta0) + p_max_post * np.cos(delta_max) - p_max_fault * np.cos(delta0)) / (p_max_post - p_max_fault)
    if np.any(cos_cr > np.cos(delta0)):
        raise ValueError("the post-fault system cannot recover even with instantaneous clearing")
    return delta0, np.arccos(np.maximum(cos_cr, np.cos(delta_max)))


def critical_clearing_time(inertia_constant, p_mech, p_max_pre, frequency=60.0, p_max_post=None):
    r"""Critical clearing time of a bolted fault that drops the electrical power to zero.

    With :math:`P_e = 0` during the fault the swing equation integrates to
    :math:`\delta = \delta_0 + \omega_s P_m t^2/4H`, so the critical
    clearing angle of :func:`equal_area_clearing_angle` is reached at

    .. math::

        t_{cr} = \sqrt{\frac{4H(\delta_{cr} - \delta_0)}{\omega_s P_m}},
        \qquad \omega_s = 2\pi f.

    Parameters
    ----------
    inertia_constant : float or array_like
        Inertia constant :math:`H` [s] (stored kinetic energy per unit
        rating), positive.
    p_mech : float or array_like
        Mechanical input power [pu], positive.
    p_max_pre : float or array_like
        Prefault power-angle amplitude [pu], above ``p_mech``.
    frequency : float or array_like, optional
        System frequency [Hz], positive.
    p_max_post : float or array_like, optional
        Post-fault amplitude [pu]; defaults to ``p_max_pre`` (the fault is
        cleared without losing a line).

    Returns
    -------
    float or ndarray
        Critical clearing time [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. For a fault that leaves some
    power flowing there is no closed form: integrate with
    :func:`swing_response`.

    Examples
    --------
    >>> round(float(critical_clearing_time(5.0, 1.0, 2.0)), 4)
    0.2142
    """
    _finite(inertia_constant=inertia_constant, frequency=frequency)
    _positive(inertia_constant=inertia_constant, frequency=frequency)
    p_max_post = p_max_pre if p_max_post is None else p_max_post
    inertia_constant, frequency, p_mech = (np.asarray(x, dtype=float) for x in (inertia_constant, frequency, p_mech))
    delta0, delta_cr = equal_area_clearing_angle(p_mech, p_max_pre, 0.0, p_max_post)
    return np.sqrt(4 * inertia_constant * (delta_cr - delta0) / (2 * np.pi * frequency * p_mech))


def swing_response(
    inertia_constant, p_mech, p_max_pre, p_max_fault, p_max_post, clearing_time, t_end, frequency=60.0, damping=0.0, n_points=2001
):
    r"""Swing curve of a generator against an infinite bus through a cleared fault.

    Integrates the swing equation

    .. math::

        \frac{2H}{\omega_s}\frac{d^2\delta}{dt^2} = P_m - P_{max}\sin\delta
        - D\frac{d\delta}{dt}

    from rest at :math:`\delta_0 = \arcsin(P_m/P_1)`, with
    :math:`P_{max} = P_2` until ``clearing_time`` and :math:`P_3` after.
    The two intervals are integrated separately, so the step change in
    :math:`P_{max}` falls on a sample.

    Parameters
    ----------
    inertia_constant : float
        Inertia constant :math:`H` [s], positive.
    p_mech : float
        Mechanical input power [pu], positive.
    p_max_pre, p_max_fault, p_max_post : float
        Power-angle amplitudes before, during and after the fault [pu]; see
        :func:`equal_area_clearing_angle`.
    clearing_time : float
        Fault duration [s], in :math:`(0, t_{end})`.
    t_end : float
        End of the simulation [s], positive.
    frequency : float, optional
        System frequency [Hz], positive.
    damping : float, optional
        Damping coefficient :math:`D` [pu power s/rad], nonnegative.
    n_points : int, optional
        Number of output samples over :math:`[0, t_{end}]`, at least three
        (the clearing instant is added to them).

    Returns
    -------
    SwingResult
        Rotor angle [rad] and speed deviation [rad/s] against time [s], with
        a first-swing ``stable`` flag.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. Integrated with
    :func:`engineeringkit.integrators.integrate_ode` (RK45). All parameters
    are scalar-only: the result is one swing curve.

    Examples
    --------
    >>> r = swing_response(5.0, 1.0, 2.0, 0.0, 2.0, 0.2, 1.5)
    >>> r.stable, round(float(np.degrees(r.delta.max())), 1)
    (True, 119.4)
    """
    _check_powers(p_mech, p_max_pre, p_max_fault, p_max_post)
    _finite(inertia_constant=inertia_constant, clearing_time=clearing_time, t_end=t_end, frequency=frequency, damping=damping)
    _positive(inertia_constant=inertia_constant, t_end=t_end, frequency=frequency)
    _nonnegative(damping=damping)
    _integer("n_points", n_points, minimum=3)
    if not 0 < clearing_time < t_end:
        raise ValueError("clearing_time must lie in (0, t_end)")
    omega_s = 2 * np.pi * frequency
    delta0 = float(np.arcsin(p_mech / p_max_pre))

    def rhs(p_max):
        return lambda _, y: [y[1], omega_s / (2 * inertia_constant) * (p_mech - p_max * np.sin(y[0]) - damping * y[1])]

    grid = np.linspace(0, t_end, n_points)
    before = np.append(grid[grid < clearing_time], clearing_time)
    after = np.insert(grid[grid > clearing_time], 0, clearing_time)
    fault = integrate_ode(rhs(p_max_fault), [delta0, 0.0], before)
    post = integrate_ode(rhs(p_max_post), fault.history[-1], after)
    t = np.concatenate((before, after[1:]))
    history = np.vstack((fault.history, post.history[1:]))
    max_stable_angle = float(np.pi - np.arcsin(p_mech / p_max_post))
    return SwingResult(t, history[:, 0], history[:, 1], float(clearing_time), max_stable_angle)
