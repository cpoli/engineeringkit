r"""Aerodynamic forces from coefficients, circulation, and thin-airfoil theory.

The force on a body in a stream is written as a dimensionless coefficient
times the dynamic pressure :math:`q = \tfrac12\rho V^2` times a reference area.
That bookkeeping goes back to Cayley and Lilienthal; the circulation theory of
Kutta and Joukowski explains *why* a wing produces lift, and Munk and
Glauert's thin-airfoil theory predicts the coefficient from section geometry.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "dynamic_pressure",
    "lift",
    "drag",
    "lift_coefficient",
    "drag_coefficient",
    "stall_speed",
    "kutta_joukowski_lift",
    "thin_airfoil_cl",
    "naca4_zero_lift_angle",
    "naca4_moment_coefficient",
    "aerodynamic_center",
]


def dynamic_pressure(rho, V):
    r"""Dynamic pressure :math:`q = \tfrac12 \rho V^2`.

    Bernoulli's kinetic-energy term per unit volume: the pressure rise
    available when a stream is brought to rest without losses.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    V : float or array_like
        Flow speed [m/s], nonnegative.

    Returns
    -------
    float or ndarray
        Dynamic pressure [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; assumptions and
    validity limits are in :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> float(dynamic_pressure(1.2, 50.0))
    1500.0
    """
    _positive(rho=rho)
    _nonnegative(V=V)
    return 0.5 * rho * np.asarray(V) ** 2


def lift(q, S, Cl):
    r"""Lift force :math:`L = q S C_L`.

    Parameters
    ----------
    q : float or array_like
        Dynamic pressure [Pa], nonnegative.
    S : float or array_like
        Reference (wing planform) area [m^2], positive.
    Cl : float or array_like
        Lift coefficient [-]; signed.

    Returns
    -------
    float or ndarray
        Lift [N], perpendicular to the free stream.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; assumptions and
    validity limits are in :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> float(lift(1500.0, 2.0, 0.5))
    1500.0
    """
    _finite(Cl=Cl)
    _nonnegative(q=q)
    _positive(S=S)
    return q * S * np.asarray(Cl)


def drag(q, S, Cd):
    r"""Drag force :math:`D = q S C_D`.

    Parameters
    ----------
    q : float or array_like
        Dynamic pressure [Pa], nonnegative.
    S : float or array_like
        Reference area [m^2], positive.
    Cd : float or array_like
        Drag coefficient [-], nonnegative.

    Returns
    -------
    float or ndarray
        Drag [N], parallel to the free stream.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; assumptions and
    validity limits are in :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> float(drag(1500.0, 2.0, 0.02))
    60.0
    """
    _nonnegative(q=q, Cd=Cd)
    _positive(S=S)
    return q * S * np.asarray(Cd)


def lift_coefficient(L, rho, V, S):
    r"""Lift coefficient :math:`C_L = L / (q S)` from a measured lift force.

    Parameters
    ----------
    L : float or array_like
        Lift [N]; signed.
    rho : float or array_like
        Density [kg/m^3], positive.
    V : float or array_like
        Speed [m/s], positive.
    S : float or array_like
        Reference area [m^2], positive.

    Returns
    -------
    float or ndarray
        Lift coefficient [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; assumptions and
    validity limits are in :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> float(lift_coefficient(1500.0, 1.2, 50.0, 2.0))
    0.5
    """
    _finite(L=L)
    _positive(rho=rho, V=V, S=S)
    return np.asarray(L) / dynamic_pressure(rho, V) / S


def drag_coefficient(D, rho, V, S):
    r"""Drag coefficient :math:`C_D = D / (q S)` from a measured drag force.

    Parameters
    ----------
    D : float or array_like
        Drag [N], nonnegative.
    rho : float or array_like
        Density [kg/m^3], positive.
    V : float or array_like
        Speed [m/s], positive.
    S : float or array_like
        Reference area [m^2], positive.

    Returns
    -------
    float or ndarray
        Drag coefficient [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; assumptions and
    validity limits are in :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> float(drag_coefficient(60.0, 1.2, 50.0, 2.0))
    0.02
    """
    _nonnegative(D=D)
    _positive(rho=rho, V=V, S=S)
    return np.asarray(D) / dynamic_pressure(rho, V) / S


def stall_speed(weight, rho, S, Clmax):
    r"""Minimum level-flight speed :math:`V_s = \sqrt{2W / (\rho S C_{L,\max})}`.

    Setting lift equal to weight at the largest attainable lift coefficient
    gives the slowest speed at which the wing can still support the aircraft.

    Parameters
    ----------
    weight : float or array_like
        Aircraft weight [N], positive.
    rho : float or array_like
        Air density [kg/m^3], positive.
    S : float or array_like
        Wing area [m^2], positive.
    Clmax : float or array_like
        Maximum lift coefficient [-], positive.

    Returns
    -------
    float or ndarray
        Stall speed [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; assumptions and
    validity limits are in :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(stall_speed(12000.0, 1.2, 20.0, 1.5)), 4)
    25.8199
    """
    _positive(weight=weight, rho=rho, S=S, Clmax=Clmax)
    return np.sqrt(2 * np.asarray(weight) / (rho * S * Clmax))


def kutta_joukowski_lift(rho, V, circulation):
    r"""Lift per unit span :math:`L' = \rho V \Gamma` (Kutta-Joukowski theorem).

    Any two-dimensional body in a steady inviscid stream carries a lift per
    unit span proportional to the circulation :math:`\Gamma` around it, and
    no drag (d'Alembert).

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3], positive.
    V : float or array_like
        Free-stream speed [m/s], nonnegative.
    circulation : float or array_like
        Circulation :math:`\Gamma` [m^2/s]; positive for clockwise flow
        about a body in a stream moving in the +x direction.

    Returns
    -------
    float or ndarray
        Lift per unit span [N/m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; assumptions and
    validity limits are in :doc:`/guides/aerodynamics`.

    References
    ----------
    M. W. Kutta, *Illustrierte Aeronautische Mitteilungen* 6 (1902);
    N. Joukowski, *Bull. Inst. Aérodyn. Koutchino* 1 (1906).

    Examples
    --------
    >>> float(kutta_joukowski_lift(1.2, 50.0, 10.0))
    600.0
    """
    _positive(rho=rho)
    _nonnegative(V=V)
    _finite(circulation=circulation)
    return rho * V * np.asarray(circulation)


def thin_airfoil_cl(alpha_rad, alpha_l0=0.0):
    r"""Thin-airfoil section lift coefficient :math:`c_l = 2\pi(\alpha - \alpha_{L=0})`.

    Parameters
    ----------
    alpha_rad : float or array_like
        Angle of attack [rad].
    alpha_l0 : float, optional
        Zero-lift angle of attack [rad]; zero for a symmetric section and
        negative for positive camber (see :func:`naca4_zero_lift_angle`).

    Returns
    -------
    float or ndarray
        Section lift coefficient [-].

    Notes
    -----
    Small angles and attached, inviscid, incompressible flow are assumed;
    stall is not predicted. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/aerodynamics`.

    References
    ----------
    H. Glauert, *The Elements of Aerofoil and Airscrew Theory* (Cambridge
    University Press, 1926), ch. 7.

    Examples
    --------
    >>> round(float(thin_airfoil_cl(0.1)), 6)
    0.628319
    """
    _finite(alpha_rad=alpha_rad, alpha_l0=alpha_l0)
    return 2 * np.pi * (np.asarray(alpha_rad) - alpha_l0)


def naca4_zero_lift_angle(m, p):
    r"""Thin-airfoil zero-lift angle of a NACA four-digit camber line.

    Glauert's integral

    .. math::

        \alpha_{L=0} = -\frac{1}{\pi}\int_0^\pi \frac{dz}{dx}(\theta)\,(\cos\theta - 1)\,d\theta,
        \qquad x = \tfrac12(1-\cos\theta),

    is evaluated numerically for the two-parabola camber line with maximum
    camber ``m`` at chordwise position ``p``.

    Parameters
    ----------
    m : float
        Maximum camber as a fraction of chord, in ``[0, 0.095]``
        (first digit / 100).
    p : float
        Position of maximum camber as a fraction of chord, in ``(0, 1)``
        (second digit / 10).

    Returns
    -------
    float
        Zero-lift angle of attack [rad]; zero for ``m = 0``.

    Notes
    -----
    ``m`` and ``p`` are scalar-only: each camber line needs its own quadrature.

    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(np.degrees(naca4_zero_lift_angle(0.02, 0.4))), 2)
    -2.08
    """
    _interval("m", m, 0, 0.095)
    _interval("p", p, 0, 1, closed=False)
    theta_p = np.arccos(1 - 2 * p)

    def front(theta):
        x = 0.5 * (1 - np.cos(theta))
        return 2 * m / p**2 * (p - x) * (np.cos(theta) - 1)

    def back(theta):
        x = 0.5 * (1 - np.cos(theta))
        return 2 * m / (1 - p) ** 2 * (p - x) * (np.cos(theta) - 1)

    integral = quad(front, 0, theta_p)[0] + quad(back, theta_p, np.pi)[0]
    return -integral / np.pi


def naca4_moment_coefficient(m, p):
    r"""Thin-airfoil quarter-chord pitching moment of a NACA four-digit camber line.

    Glauert's Fourier coefficients of the camber-line slope,

    .. math::

        A_n = \frac{2}{\pi}\int_0^\pi \frac{dz}{dx}(\theta)\cos n\theta\,d\theta,

    give :math:`c_{m,c/4} = \tfrac{\pi}{4}(A_2 - A_1)`, which does not depend
    on the angle of attack: the quarter chord is the thin-airfoil
    aerodynamic centre.

    Parameters
    ----------
    m : float
        Maximum camber as a fraction of chord, in ``[0, 0.095]``
        (first digit / 100).
    p : float
        Position of maximum camber as a fraction of chord, in ``(0, 1)``
        (second digit / 10).

    Returns
    -------
    float
        Moment coefficient about the quarter chord [-], nose up positive;
        negative for positive camber and zero for a symmetric section.

    Notes
    -----
    The moment about the leading edge follows as
    :math:`c_{m,le} = c_{m,c/4} - c_l/4`. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Fundamentals of Aerodynamics*, 6th ed. (McGraw-Hill,
    2017), sec. 4.8 and example 4.6.

    Examples
    --------
    >>> round(float(naca4_moment_coefficient(0.02, 0.4)), 4)
    -0.0531
    """
    _interval("m", m, 0, 0.095)
    _interval("p", p, 0, 1, closed=False)
    theta_p = np.arccos(1 - 2 * p)

    def coefficient(n):
        def front(theta):
            return 2 * m / p**2 * (p - 0.5 * (1 - np.cos(theta))) * np.cos(n * theta)

        def back(theta):
            return 2 * m / (1 - p) ** 2 * (p - 0.5 * (1 - np.cos(theta))) * np.cos(n * theta)

        return 2 / np.pi * (quad(front, 0, theta_p)[0] + quad(back, theta_p, np.pi)[0])

    return np.pi / 4 * (coefficient(2) - coefficient(1))


def aerodynamic_center(a0, m0):
    r"""Chordwise position of the aerodynamic centre, :math:`x_{ac}/c = \tfrac14 - m_0/a_0`.

    The aerodynamic centre is the point about which the pitching moment does
    not change with angle of attack.

    Parameters
    ----------
    a0 : float or array_like
        Section lift-curve slope :math:`dc_l/d\alpha` [1/rad], positive.
    m0 : float or array_like
        Slope of the quarter-chord moment coefficient,
        :math:`dc_{m,c/4}/d\alpha` [1/rad]; zero in thin-airfoil theory.

    Returns
    -------
    float or ndarray
        Aerodynamic-centre position as a fraction of chord from the leading
        edge [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Fundamentals of Aerodynamics*, 6th ed. (McGraw-Hill,
    2017), sec. 4.9.

    Examples
    --------
    >>> round(float(aerodynamic_center(6.0, -0.06)), 4)
    0.26
    """
    _positive(a0=a0)
    _finite(m0=m0)
    return 0.25 - np.asarray(m0) / a0
