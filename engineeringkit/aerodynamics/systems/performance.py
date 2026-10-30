r"""Airplane performance: climb, level turns, take-off and propeller range.

Each result is a force or energy balance on a point-mass aircraft. Excess
thrust lifts it (climb), extra lift curves its path (turn), and net thrust
over the ground roll accelerates it to lift-off. Burning fuel at a rate set
by shaft power gives Breguet's propeller range.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.aerodynamics.systems.wings import STANDARD_GRAVITY

__all__ = [
    "rate_of_climb",
    "climb_angle",
    "load_factor",
    "turn_radius",
    "turn_rate",
    "takeoff_ground_roll",
    "breguet_range_propeller",
]


def rate_of_climb(thrust, drag, V, weight):
    r"""Steady rate of climb from excess power, :math:`R/C = (T - D)V/W`.

    Parameters
    ----------
    thrust : float or array_like
        Thrust available [N], nonnegative.
    drag : float or array_like
        Drag at the flight condition [N], nonnegative.
    V : float or array_like
        True airspeed [m/s], nonnegative.
    weight : float or array_like
        Aircraft weight [N], positive.

    Returns
    -------
    float or ndarray
        Vertical speed [m/s]; negative when drag exceeds thrust.

    Notes
    -----
    Uses the drag of level flight, a good approximation for climb angles
    below about 15°. Inputs are raw SI numbers with no dimensional metadata;
    see :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Introduction to Flight*, 8th ed. (McGraw-Hill, 2016), sec. 6.8.

    Examples
    --------
    >>> float(rate_of_climb(3000.0, 1500.0, 60.0, 15000.0))
    6.0
    """
    _nonnegative(thrust=thrust, drag=drag, V=V)
    _positive(weight=weight)
    return (np.asarray(thrust) - drag) * V / weight


def climb_angle(thrust, drag, weight):
    r"""Steady climb angle :math:`\gamma = \arcsin\left[(T - D)/W\right]`.

    Parameters
    ----------
    thrust : float or array_like
        Thrust available [N], nonnegative.
    drag : float or array_like
        Drag [N], nonnegative.
    weight : float or array_like
        Aircraft weight [N], positive, with :math:`|T - D| \le W`.

    Returns
    -------
    float or ndarray
        Flight-path angle [rad].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(np.degrees(climb_angle(3000.0, 1500.0, 15000.0))), 3)
    5.739
    """
    _nonnegative(thrust=thrust, drag=drag)
    _positive(weight=weight)
    ratio = (np.asarray(thrust) - drag) / weight
    _interval("(thrust - drag) / weight", ratio, -1, 1)
    return np.arcsin(ratio)


def load_factor(bank_angle):
    r"""Load factor of a level, coordinated turn, :math:`n = L/W = 1/\cos\phi`.

    Parameters
    ----------
    bank_angle : float or array_like
        Bank angle :math:`\phi` [rad], in ``[0, pi/2)``.

    Returns
    -------
    float or ndarray
        Load factor [-]; 2 at a 60° bank.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(load_factor(np.radians(60.0))), 6)
    2.0
    """
    _interval("bank_angle", bank_angle, 0, np.pi / 2, closed="left")
    return 1 / np.cos(bank_angle)


def _turn(V, n, g):
    _nonnegative(V=V)
    _interval("n", n, 1, np.inf, closed=False)
    _positive(g=g)
    return np.sqrt(np.asarray(n) ** 2 - 1)


def turn_radius(V, n, g=STANDARD_GRAVITY):
    r"""Radius of a level turn, :math:`R = V^2/(g\sqrt{n^2 - 1})`.

    Parameters
    ----------
    V : float or array_like
        True airspeed [m/s], nonnegative.
    n : float or array_like
        Load factor [-], greater than one.
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    float or ndarray
        Turn radius [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Introduction to Flight*, 8th ed. (McGraw-Hill, 2016), sec. 6.17.

    Examples
    --------
    >>> round(float(turn_radius(60.0, 2.0)), 2)
    211.94
    """
    return np.asarray(V) ** 2 / (g * _turn(V, n, g))


def turn_rate(V, n, g=STANDARD_GRAVITY):
    r"""Turn rate of a level turn, :math:`\omega = g\sqrt{n^2 - 1}/V`.

    Parameters
    ----------
    V : float or array_like
        True airspeed [m/s], positive.
    n : float or array_like
        Load factor [-], greater than one.
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    float or ndarray
        Turn rate [rad/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(np.degrees(turn_rate(60.0, 2.0))), 3)
    16.22
    """
    _positive(V=V)
    return g * _turn(V, n, g) / np.asarray(V)


def takeoff_ground_roll(weight, rho, S, Clmax, thrust, resistance=0.0, g=STANDARD_GRAVITY):
    r"""Take-off ground roll to lift-off at :math:`1.2V_{stall}`.

    .. math::

        s_{LO} = \frac{1.44\,W^2}{g\,\rho\,S\,C_{L,\max}\,(T - R)}

    with the net accelerating force :math:`T - R` taken constant along the
    roll, where :math:`R = D + \mu_r(W - L)` is the drag and rolling friction
    evaluated at :math:`0.7V_{LO}`.

    Parameters
    ----------
    weight : float or array_like
        Take-off weight [N], positive.
    rho : float or array_like
        Air density [kg/m^3], positive.
    S : float or array_like
        Wing area [m^2], positive.
    Clmax : float or array_like
        Maximum lift coefficient in take-off configuration [-], positive.
    thrust : float or array_like
        Mean thrust over the roll [N], positive.
    resistance : float or array_like, optional
        Mean drag plus rolling friction [N], nonnegative and less than ``thrust``.
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    float or ndarray
        Ground-roll distance [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Introduction to Flight*, 8th ed. (McGraw-Hill, 2016), sec. 6.15.

    Examples
    --------
    >>> round(float(takeoff_ground_roll(10.7e3, 1.225, 16.2, 1.6, 2400.0, 400.0)), 1)
    264.7
    """
    _positive(weight=weight, rho=rho, S=S, Clmax=Clmax, thrust=thrust, g=g)
    _nonnegative(resistance=resistance)
    net = np.asarray(thrust) - resistance
    if np.any(net <= 0):
        raise ValueError("thrust must exceed resistance")
    return 1.44 * np.asarray(weight) ** 2 / (g * rho * S * Clmax * net)


def breguet_range_propeller(eta_prop, psfc, lift_to_drag, W_initial, W_final, g=STANDARD_GRAVITY):
    r"""Propeller-aircraft cruise range :math:`R = \dfrac{\eta_{pr}}{g\,c_P}\dfrac{L}{D}\ln\dfrac{W_i}{W_f}` (Breguet).

    Shaft power :math:`P = DV/\eta_{pr}` burns fuel at :math:`\dot m = c_P P`
    and drag is :math:`W/(L/D)`, so the speed cancels: a propeller
    aircraft's range depends on speed only through :math:`L/D`.

    Parameters
    ----------
    eta_prop : float or array_like
        Propeller efficiency [-], in ``(0, 1]``.
    psfc : float or array_like
        Power-specific fuel consumption [kg/(W s)], positive
        (``7.6e-8``, about 0.45 lb/(hp h), for a piston engine).
    lift_to_drag : float or array_like
        Cruise lift-to-drag ratio [-], positive.
    W_initial, W_final : float or array_like
        Weight at the start and end of cruise [N], with ``W_initial >= W_final > 0``.
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    float or ndarray
        Still-air range [m]; largest at the maximum-:math:`L/D` speed.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Introduction to Flight*, 8th ed. (McGraw-Hill, 2016),
    sec. 6.12.

    Examples
    --------
    >>> round(float(breguet_range_propeller(0.8, 7.6e-8, 13.6, 13.1e3, 11.7e3)) / 1e3)
    1650
    """
    _interval("eta_prop", eta_prop, 0, 1, closed="right")
    _positive(psfc=psfc, lift_to_drag=lift_to_drag, W_final=W_final, g=g)
    _finite(W_initial=W_initial)
    if np.any(np.asarray(W_initial) < W_final):
        raise ValueError("W_initial must be at least W_final")
    return np.asarray(eta_prop) / (g * psfc) * lift_to_drag * np.log(np.asarray(W_initial) / W_final)
