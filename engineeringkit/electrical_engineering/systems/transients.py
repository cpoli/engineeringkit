r"""First- and second-order transients: RC and RL responses and the series RLC step.

Helmholtz showed that the current in an inductive circuit rises to its
final value on the time scale :math:`L/R`, and Kelvin's analysis of the
first submarine telegraph cables showed that a capacitance charged through
a resistance responds on the time scale :math:`RC`. With both an inductance
and a capacitance the circuit is second order, and Kelvin found that a
Leyden jar discharges by oscillation when :math:`R < 2\sqrt{L/C}`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["rc_time_constant", "rl_time_constant", "rc_charging_voltage", "rl_current", "rlc_damping", "rlc_step_voltage"]


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


def rl_current(v_source, R, L, t, i0=0.0):
    r"""Inductor current :math:`i(t) = V_s/R + (i_0 - V_s/R)e^{-Rt/L}` after a step.

    Parameters
    ----------
    v_source : float
        Constant source voltage applied at :math:`t = 0` [V].
    R : float
        Series resistance [ohm], positive.
    L : float
        Inductance [H], positive.
    t : float or array_like
        Time since the step [s], nonnegative.
    i0 : float, optional
        Initial inductor current [A].

    Returns
    -------
    ndarray
        Inductor current [A], shaped like ``t``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    H. Helmholtz, "Ueber die Dauer und den Verlauf der durch
    Stromesschwankungen inducirten elektrischen Ströme," *Ann. Phys.* 83
    (1851) 505-540.

    Examples
    --------
    >>> round(float(rl_current(10.0, 5.0, 0.5, 0.1)), 4)
    1.2642
    """
    _finite(v_source=v_source, R=R, L=L, t=t, i0=i0)
    _nonnegative(t=t)
    t = np.asarray(t, float)
    final = v_source / R
    return final + (i0 - final) * np.exp(-t / rl_time_constant(L, R))


def rlc_damping(R, L, C):
    r"""Neper frequency, resonant frequency and damping regime of a series RLC circuit.

    The natural response of a series RLC circuit has the characteristic
    roots :math:`s = -\alpha \pm \sqrt{\alpha^2 - \omega_0^2}` with
    :math:`\alpha = R/2L` and :math:`\omega_0 = 1/\sqrt{LC}`. It is
    overdamped for :math:`\alpha > \omega_0`, critically damped for
    :math:`\alpha = \omega_0` (:math:`R = 2\sqrt{L/C}`) and underdamped,
    an oscillation decaying as :math:`e^{-\alpha t}`, for
    :math:`\alpha < \omega_0`.

    Parameters
    ----------
    R : float
        Series resistance [ohm], nonnegative.
    L : float
        Inductance [H], positive.
    C : float
        Capacitance [F], positive.

    Returns
    -------
    tuple of (float, float, str)
        ``(alpha, omega0, regime)``: the neper frequency [1/s], the
        undamped resonant frequency [rad/s] and one of ``"overdamped"``,
        ``"critically damped"`` (within a relative 1e-9), ``"underdamped"``
        or ``"undamped"`` (:math:`R = 0`).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. The damping ratio is
    :math:`\zeta = \alpha/\omega_0`. ``R``, ``L`` and ``C`` are
    scalar-only: the regime names one circuit.

    References
    ----------
    W. Thomson, "On transient electric currents," *Phil. Mag.* 5 (1853)
    393-405.

    Examples
    --------
    >>> rlc_damping(2.0, 1.0, 1.0)
    (1.0, 1.0, 'critically damped')
    """
    _finite(R=R, L=L, C=C)
    _nonnegative(R=R)
    _positive(L=L, C=C)
    alpha, omega0 = R / (2 * L), 1 / np.sqrt(L * C)
    if alpha == 0:
        regime = "undamped"
    elif np.isclose(alpha, omega0, rtol=1e-9, atol=0):
        regime = "critically damped"
    else:
        regime = "overdamped" if alpha > omega0 else "underdamped"
    return float(alpha), float(omega0), regime


def rlc_step_voltage(v_source, R, L, C, t):
    r"""Capacitor voltage of a series RLC circuit after a voltage step, from rest.

    With :math:`\alpha` and :math:`\omega_0` from :func:`rlc_damping`,

    .. math::

        \frac{v_C(t)}{V_s} = \begin{cases}
        1 + \dfrac{s_2 e^{s_1 t} - s_1 e^{s_2 t}}{s_1 - s_2} & \text{overdamped},\\
        1 - (1 + \alpha t)e^{-\alpha t} & \text{critically damped},\\
        1 - e^{-\alpha t}\left(\cos\omega_d t + \dfrac{\alpha}{\omega_d}
        \sin\omega_d t\right) & \text{underdamped},
        \end{cases}

    with :math:`s_{1,2} = -\alpha \pm \sqrt{\alpha^2 - \omega_0^2}` and
    :math:`\omega_d = \sqrt{\omega_0^2 - \alpha^2}`.

    Parameters
    ----------
    v_source : float
        Constant source voltage applied at :math:`t = 0` [V].
    R : float
        Series resistance [ohm], nonnegative.
    L : float
        Inductance [H], positive.
    C : float
        Capacitance [F], positive.
    t : float or array_like
        Time since the step [s], nonnegative.

    Returns
    -------
    ndarray
        Capacitor voltage [V], shaped like ``t``. The loop current is
        :math:`C\,dv_C/dt`, as returned by
        :meth:`SeriesRLC.step_response <engineeringkit.electrical_engineering.SeriesRLC.step_response>`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. Within the relative 1e-9 band
    that :func:`rlc_damping` calls critical, the critical formula is used.
    ``R``, ``L`` and ``C`` are scalar-only (one circuit, swept over ``t``).

    Examples
    --------
    >>> round(float(rlc_step_voltage(1.0, 2.0, 1.0, 1.0, 1.0)), 4)
    0.2642
    """
    _finite(v_source=v_source, t=t)
    _nonnegative(t=t)
    alpha, omega0, regime = rlc_damping(R, L, C)
    t = np.asarray(t, float)
    if regime == "overdamped":
        root = np.sqrt(alpha**2 - omega0**2)
        s1, s2 = -alpha + root, -alpha - root
        shape = 1 + (s2 * np.exp(s1 * t) - s1 * np.exp(s2 * t)) / (s1 - s2)
    elif regime == "critically damped":
        shape = 1 - (1 + alpha * t) * np.exp(-alpha * t)
    else:
        wd = np.sqrt(omega0**2 - alpha**2)
        shape = 1 - np.exp(-alpha * t) * (np.cos(wd * t) + alpha / wd * np.sin(wd * t))
    return v_source * shape
