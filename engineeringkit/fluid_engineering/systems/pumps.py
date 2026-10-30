r"""Pumps: Euler's head, pump and system curves, the affinity laws and cavitation.

Euler (1754) found the head a rotor gives the fluid from the change in its
angular momentum, :math:`H = (u_2 V_{t2} - u_1 V_{t1})/g`. A real pump's
head falls with flow along its curve, the system needs a head that rises
with flow, and the pump runs where the two cross. The affinity laws scale
a curve to another speed or impeller diameter. A pump adds head :math:`H` to a flow :math:`Q`, so the power given to the
fluid is :math:`\rho g Q H`. Where the local pressure falls to the vapour
pressure the liquid boils and the collapsing bubbles pit the blades. Thoma
(1924) found that turbines of one design all start to cavitate at the same
ratio :math:`\sigma = \mathrm{NPSH}/H` of the suction head margin to the
machine head.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.core.base import OperatingPoint
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G

__all__ = [
    "pump_hydraulic_power",
    "pump_input_power",
    "pump_specific_speed",
    "cavitation_number",
    "npsh_available",
    "thoma_sigma",
    "euler_turbomachine_head",
    "pump_curve_head",
    "system_curve_head",
    "pump_operating_point",
    "affinity_laws",
]


def pump_hydraulic_power(rho, Q, head, g=G):
    r"""Power delivered to the fluid, :math:`P = \rho g Q H`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    Q : float or array_like
        Volume flow [m^3/s], nonnegative.
    head : float or array_like
        Head added by the pump [m], nonnegative.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Hydraulic power [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> pump_hydraulic_power(1000.0, 0.01, 10.0, g=10.0)
    1000.0
    """
    _finite(rho=rho, Q=Q, head=head, g=g)
    _positive(rho=rho, g=g)
    _nonnegative(Q=Q, head=head)
    return rho * g * Q * head


def pump_input_power(rho, Q, head, efficiency, g=G):
    r"""Shaft power a pump needs, :math:`P = \rho g Q H/\eta`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    Q : float or array_like
        Volume flow [m^3/s], nonnegative.
    head : float or array_like
        Head added by the pump [m], nonnegative.
    efficiency : float or array_like
        Pump efficiency [-], in (0, 1].
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Shaft power [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> pump_input_power(1000.0, 0.01, 10.0, 0.5, g=10.0)
    2000.0
    """
    _finite(rho=rho, Q=Q, head=head, efficiency=efficiency, g=g)
    _interval("efficiency", efficiency, 0, 1)
    if np.any(np.asarray(efficiency) == 0):
        raise ValueError("efficiency must be positive")
    return pump_hydraulic_power(rho, Q, head, g) / efficiency


def pump_specific_speed(omega, Q, head, g=G):
    r"""Dimensionless specific speed :math:`\Omega_s = \omega\sqrt Q/(gH)^{3/4}`.

    Parameters
    ----------
    omega : float or array_like
        Shaft speed [rad/s], nonnegative.
    Q : float or array_like
        Volume flow at best efficiency [m^3/s], nonnegative.
    head : float or array_like
        Head at best efficiency [m], positive.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Specific speed [-]: below about 1 for radial pumps, above about 4
        for axial ones.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Industry values in rpm, gpm and feet
    differ by a constant factor.

    Examples
    --------
    >>> float(pump_specific_speed(10.0, 4.0, 1.0, g=1.0))
    20.0
    """
    _finite(omega=omega, Q=Q, head=head, g=g)
    _positive(head=head, g=g)
    _nonnegative(Q=Q, omega=omega)
    return omega * np.sqrt(Q) / (g * head) ** 0.75


def cavitation_number(p, p_vapor, rho, v):
    r"""Cavitation number :math:`\mathrm{Ca} = (p - p_v)/(\tfrac12\rho v^2)`.

    Parameters
    ----------
    p : float or array_like
        Reference absolute pressure [Pa].
    p_vapor : float or array_like
        Vapour pressure of the liquid [Pa].
    rho : float or array_like
        Liquid density [kg/m^3], positive.
    v : float or array_like
        Reference speed [m/s], nonzero.

    Returns
    -------
    float or ndarray
        Cavitation number [-]; negative when the pressure is below the
        vapour pressure.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> cavitation_number(3000.0, 1000.0, 1000.0, 2.0)
    1.0
    """
    _finite(p=p, p_vapor=p_vapor, rho=rho, v=v)
    _positive(rho=rho)
    if np.any(np.asarray(v) == 0):
        raise ValueError("velocity must be nonzero")
    return (p - p_vapor) / (0.5 * rho * v * v)


def npsh_available(p_inlet, v_inlet, p_vapor, rho, g=G):
    r"""Net positive suction head available at a pump inlet.

    .. math::

        \mathrm{NPSH}_A = \frac{p - p_v}{\rho g} + \frac{v^2}{2g}

    Parameters
    ----------
    p_inlet : float or array_like
        Absolute static pressure at the pump inlet [Pa], positive.
    v_inlet : float or array_like
        Mean velocity at the inlet [m/s].
    p_vapor : float or array_like
        Vapour pressure of the liquid [Pa], nonnegative.
    rho : float or array_like
        Liquid density [kg/m^3], positive.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        NPSH available [m]; the pump cavitates when it falls below the
        NPSH it requires.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. For a pump drawing from an open sump,
    :math:`p/(\rho g) + v^2/(2g) = p_{atm}/(\rho g) - z_s - h_{L}`, the
    atmosphere less the suction lift and the suction-line losses.

    Examples
    --------
    >>> npsh_available(103000.0, 0.0, 3000.0, 1000.0, g=10.0)
    10.0
    """
    _finite(p_inlet=p_inlet, v_inlet=v_inlet, p_vapor=p_vapor, rho=rho, g=g)
    _positive(p_inlet=p_inlet, rho=rho, g=g)
    _nonnegative(p_vapor=p_vapor)
    return (p_inlet - p_vapor) / (rho * g) + v_inlet * v_inlet / (2 * g)


def thoma_sigma(npsh, head):
    r"""Thoma's cavitation parameter :math:`\sigma = \mathrm{NPSH}/H`.

    Parameters
    ----------
    npsh : float or array_like
        Net positive suction head [m].
    head : float or array_like
        Head of the machine [m], positive.

    Returns
    -------
    float or ndarray
        Thoma's :math:`\sigma` [-]. A machine runs free of cavitation while
        its plant :math:`\sigma` stays above the critical value found by
        test for its design.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    References
    ----------
    D. Thoma, "Die Kavitation bei Wasserturbinen," *Wasserkraft-Jahrbuch*
    (Munich, 1924).

    Examples
    --------
    >>> thoma_sigma(5.0, 50.0)
    0.1
    """
    _finite(npsh=npsh, head=head)
    _positive(head=head)
    return npsh / head


def euler_turbomachine_head(u2, v_t2, u1=0.0, v_t1=0.0, g=G):
    r"""Ideal head of a rotor from Euler's turbomachine equation.

    .. math::

        H = \frac{u_2 V_{t2} - u_1 V_{t1}}{g}

    Parameters
    ----------
    u2 : float or array_like
        Blade speed at the outlet, :math:`\omega r_2` [m/s].
    v_t2 : float or array_like
        Tangential (swirl) component of the absolute velocity at the
        outlet [m/s].
    u1 : float or array_like, optional
        Blade speed at the inlet [m/s].
    v_t1 : float or array_like, optional
        Swirl velocity at the inlet [m/s]; zero for a pump without inlet
        guide vanes.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Ideal head [m]: positive for a pump, negative for a turbine (which
        takes head from the fluid).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. For a radial impeller with outlet
    blade angle :math:`\beta_2` and width :math:`b_2`,
    :math:`V_{t2} = u_2 - Q\cot\beta_2/(2\pi r_2 b_2)`, so the ideal head
    falls linearly with flow for backward-curved blades. Real pumps lose
    more to slip, friction and shock.

    References
    ----------
    L. Euler, "Théorie plus complette des machines qui sont mises en
    mouvement par la réaction de l'eau," *Mém. Acad. Sci. Berlin* 10
    (1754, publ. 1756) 227-295 (Eneström E222).

    Examples
    --------
    >>> euler_turbomachine_head(20.0, 10.0, g=10.0)
    20.0
    """
    _finite(u2=u2, v_t2=v_t2, u1=u1, v_t1=v_t1, g=g)
    _positive(g=g)
    return (u2 * v_t2 - u1 * v_t1) / g


def pump_curve_head(Q, shutoff_head, a, b=0.0):
    r"""Head of a pump on a quadratic curve, :math:`H = H_0 + bQ - aQ^2`.

    Parameters
    ----------
    Q : float or array_like
        Volume flow [m^3/s], nonnegative.
    shutoff_head : float or array_like
        Head at zero flow :math:`H_0` [m], positive.
    a : float or array_like
        Quadratic coefficient [s^2/m^5], nonnegative.
    b : float or array_like, optional
        Linear coefficient [s/m^2]; a small positive value gives the
        rising-then-falling curve of some radial pumps.

    Returns
    -------
    float or ndarray
        Pump head [m]; negative beyond the flow at which the pump runs out.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Fit ``shutoff_head``, ``a`` and
    ``b`` to three points of a manufacturer's curve.

    Examples
    --------
    >>> pump_curve_head(0.1, 40.0, 1000.0)
    30.0
    """
    _finite(Q=Q, shutoff_head=shutoff_head, a=a, b=b)
    _nonnegative(Q=Q, a=a)
    _positive(shutoff_head=shutoff_head)
    return shutoff_head + b * Q - a * Q * Q


def system_curve_head(Q, static_head, k):
    r"""Head a piping system needs, :math:`H = H_s + kQ^2`.

    Parameters
    ----------
    Q : float or array_like
        Volume flow [m^3/s], nonnegative.
    static_head : float or array_like
        Static lift, outlet minus inlet surface elevation plus any
        pressure difference as head [m]; may be negative.
    k : float or array_like
        Loss coefficient [s^2/m^5], nonnegative: the sum of
        :func:`pipe_resistance` over the pipes in series plus the minor
        losses as :math:`K/(2gA^2)`.

    Returns
    -------
    float or ndarray
        Required head [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Holding ``k`` fixed treats the
    friction factor as constant, which is exact for fully rough flow.

    Examples
    --------
    >>> system_curve_head(0.1, 10.0, 1000.0)
    20.0
    """
    _finite(Q=Q, static_head=static_head, k=k)
    _nonnegative(Q=Q, k=k)
    return static_head + k * Q * Q


def pump_operating_point(shutoff_head, a, static_head, k, b=0.0, series=1, parallel=1):
    r"""Flow and head where identical pumps meet a system curve.

    With :math:`n_s` pumps in series and :math:`n_p` such strings in
    parallel, the combined curve is
    :math:`H = n_s\,(H_0 + bQ/n_p - a(Q/n_p)^2)`. Setting it equal to
    :math:`H_s + kQ^2` gives a quadratic in :math:`Q` with one positive
    root.

    Parameters
    ----------
    shutoff_head : float
        Shut-off head of one pump [m], positive.
    a : float
        Quadratic coefficient of one pump [s^2/m^5], nonnegative.
    static_head : float
        Static lift of the system [m]; must be below the combined
        shut-off head.
    k : float
        System loss coefficient [s^2/m^5], nonnegative; ``a`` and ``k``
        must not both be zero.
    b : float, optional
        Linear coefficient of one pump [s/m^2].
    series : int, optional
        Pumps in series, at least 1.
    parallel : int, optional
        Strings in parallel, at least 1.

    Returns
    -------
    OperatingPoint
        Flow [m^3/s] and head [m] at the crossing, with both curves.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Series pumps add head at a given
    flow; parallel pumps add flow at a given head. Which gains more depends
    on how steep the system curve is.

    Examples
    --------
    >>> pump_operating_point(40.0, 1000.0, 10.0, 2000.0).summary()["flow_m3_s"]
    0.1
    """
    _finite(shutoff_head=shutoff_head, a=a, static_head=static_head, k=k, b=b)
    _positive(shutoff_head=shutoff_head)
    _nonnegative(a=a, k=k)
    _integer("series", series, 1)
    _integer("parallel", parallel, 1)
    A = series * a / parallel**2 + k
    B = series * b / parallel
    C = static_head - series * shutoff_head
    if A == 0:
        raise ValueError("a and k must not both be zero")
    if C >= 0:
        raise ValueError("static_head must be below the combined shut-off head")
    Q = (B + np.sqrt(B * B - 4 * A * C)) / (2 * A)
    return OperatingPoint(
        flow=float(Q),
        head=float(static_head + k * Q * Q),
        shutoff_head=float(shutoff_head),
        a=float(a),
        b=float(b),
        static_head=float(static_head),
        k=float(k),
        series=int(series),
        parallel=int(parallel),
    )


def affinity_laws(Q, head, power, speed_ratio, diameter_ratio=1.0):
    r"""Scale a pump's flow, head and power to a new speed or impeller size.

    .. math::

        \frac{Q_2}{Q_1} = \frac{N_2}{N_1}\left(\frac{D_2}{D_1}\right)^3,\quad
        \frac{H_2}{H_1} = \left(\frac{N_2}{N_1}\right)^2\left(\frac{D_2}{D_1}\right)^2,\quad
        \frac{P_2}{P_1} = \left(\frac{N_2}{N_1}\right)^3\left(\frac{D_2}{D_1}\right)^5

    Parameters
    ----------
    Q : float or array_like
        Volume flow at the original condition [m^3/s], nonnegative.
    head : float or array_like
        Head at the original condition [m].
    power : float or array_like
        Shaft power at the original condition [W], nonnegative.
    speed_ratio : float or array_like
        :math:`N_2/N_1` [-], positive.
    diameter_ratio : float or array_like, optional
        :math:`D_2/D_1` [-], positive.

    Returns
    -------
    tuple of float or ndarray
        ``(Q2, head2, power2)`` in m^3/s, m and W.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The laws hold between geometrically
    similar pumps at equal flow coefficient :math:`Q/(ND^3)`, so a point
    on one curve maps to a point of equal efficiency on the other.
    Trimming an impeller is not exact geometric similarity.

    Examples
    --------
    >>> affinity_laws(0.1, 30.0, 40e3, 2.0)
    (0.2, 120.0, 320000.0)
    """
    _finite(Q=Q, head=head, power=power, speed_ratio=speed_ratio, diameter_ratio=diameter_ratio)
    _nonnegative(Q=Q, power=power)
    _positive(speed_ratio=speed_ratio, diameter_ratio=diameter_ratio)
    n, d = speed_ratio, diameter_ratio
    return Q * n * d**3, head * n**2 * d**2, power * n**3 * d**5
