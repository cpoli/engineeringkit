r"""Finite wings: lifting-line slope, induced drag, drag polars, and cruise performance.

Prandtl's lifting-line theory replaces a wing by a bound vortex whose
trailing vortices induce a downwash; the downwash tilts the local lift
backwards, producing the induced drag :math:`C_L^2/(\pi e AR)` and lowering
the lift-curve slope relative to the two-dimensional section.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "aspect_ratio",
    "induced_drag_coefficient",
    "drag_polar",
    "finite_wing_lift_slope",
    "finite_wing_cl",
    "oswald_efficiency_estimate",
    "glide_ratio",
    "best_glide_cl",
    "breguet_range",
]

STANDARD_GRAVITY = 9.80665


def _span_efficiency(e):
    _interval("e", e, 0, 1, closed="right")


def aspect_ratio(span, area):
    r"""Wing aspect ratio :math:`AR = b^2 / S`.

    Parameters
    ----------
    span : float or array_like
        Wing span [m], positive.
    area : float or array_like
        Planform area [m^2], positive.

    Returns
    -------
    float or ndarray
        Aspect ratio [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> float(aspect_ratio(10.0, 20.0))
    5.0
    """
    _positive(span=span, area=area)
    return np.asarray(span) ** 2 / area


def induced_drag_coefficient(Cl, AR, e=1.0):
    r"""Induced (vortex) drag coefficient :math:`C_{D,i} = C_L^2/(\pi e AR)`.

    Parameters
    ----------
    Cl : float or array_like
        Wing lift coefficient [-].
    AR : float or array_like
        Aspect ratio [-], positive.
    e : float or array_like, optional
        Span efficiency factor [-] in ``(0, 1]``; one for an elliptic load.

    Returns
    -------
    float or ndarray
        Induced drag coefficient [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    L. Prandtl, "Tragflügeltheorie," *Nachr. Ges. Wiss. Göttingen* (1918).

    Examples
    --------
    >>> round(float(induced_drag_coefficient(0.5, 8.0)), 6)
    0.009947
    """
    _finite(Cl=Cl)
    _positive(AR=AR)
    _span_efficiency(e)
    return np.asarray(Cl) ** 2 / (np.pi * e * AR)


def drag_polar(Cl, Cd0, AR, e=1.0):
    r"""Parabolic drag polar :math:`C_D = C_{D0} + C_L^2/(\pi e AR)`.

    Parameters
    ----------
    Cl : float or array_like
        Lift coefficient [-].
    Cd0 : float or array_like
        Zero-lift (parasite) drag coefficient [-], nonnegative.
    AR : float or array_like
        Aspect ratio [-], positive.
    e : float or array_like, optional
        Oswald efficiency factor [-] in ``(0, 1]``.

    Returns
    -------
    float or ndarray
        Total drag coefficient [-].

    Notes
    -----
    Compressibility and Reynolds-number variation of :math:`C_{D0}` are not
    modelled. Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(drag_polar(0.5, 0.02, 8.0, 0.8)), 6)
    0.032434
    """
    _nonnegative(Cd0=Cd0)
    return Cd0 + induced_drag_coefficient(Cl, AR, e)


def finite_wing_lift_slope(AR, e=1.0, a0=2 * np.pi):
    r"""Lifting-line lift-curve slope :math:`a = a_0 / (1 + a_0/(\pi e AR))`.

    Parameters
    ----------
    AR : float or array_like
        Aspect ratio [-], positive.
    e : float or array_like, optional
        Span efficiency factor [-] in ``(0, 1]``.
    a0 : float, optional
        Two-dimensional section lift slope [1/rad]; :math:`2\pi` from thin-airfoil theory.

    Returns
    -------
    float or ndarray
        Wing lift-curve slope [1/rad]; tends to ``a0`` as ``AR`` grows.

    Notes
    -----
    Valid for unswept, moderate-to-high aspect-ratio wings in attached
    incompressible flow. Inputs are raw SI numbers; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(finite_wing_lift_slope(8.0)), 4)
    5.0265
    """
    _positive(AR=AR, a0=a0)
    _span_efficiency(e)
    return a0 / (1 + a0 / (np.pi * e * np.asarray(AR)))


def finite_wing_cl(alpha_rad, AR, e=1.0, alpha_l0=0.0):
    r"""Wing lift coefficient :math:`C_L = a(\alpha - \alpha_{L=0})` with the lifting-line slope.

    Parameters
    ----------
    alpha_rad : float or array_like
        Angle of attack [rad].
    AR : float or array_like
        Aspect ratio [-], positive.
    e : float, optional
        Span efficiency factor [-] in ``(0, 1]``.
    alpha_l0 : float, optional
        Zero-lift angle [rad].

    Returns
    -------
    float or ndarray
        Wing lift coefficient [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(finite_wing_cl(0.1, 8.0)), 5)
    0.50265
    """
    _finite(alpha_rad=alpha_rad, alpha_l0=alpha_l0)
    return finite_wing_lift_slope(AR, e) * (np.asarray(alpha_rad) - alpha_l0)


def oswald_efficiency_estimate(AR):
    r"""Empirical Oswald factor for straight wings, :math:`e = 1.78(1 - 0.045 AR^{0.68}) - 0.64`.

    Parameters
    ----------
    AR : float or array_like
        Aspect ratio [-] in ``[2, 10]``, the range of the underlying data.

    Returns
    -------
    float or ndarray
        Oswald efficiency factor [-].

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    References
    ----------
    D. P. Raymer, *Aircraft Design: A Conceptual Approach*, 6th ed. (AIAA,
    2018), eq. 12.48.

    Examples
    --------
    >>> round(float(oswald_efficiency_estimate(8.0)), 4)
    0.8106
    """
    _interval("AR", AR, 2, 10)
    return 1.78 * (1 - 0.045 * np.asarray(AR) ** 0.68) - 0.64


def glide_ratio(Cl, Cd):
    r"""Lift-to-drag ratio :math:`L/D = C_L/C_D`, the still-air glide distance per unit height lost.

    Parameters
    ----------
    Cl : float or array_like
        Lift coefficient [-].
    Cd : float or array_like
        Drag coefficient [-], positive.

    Returns
    -------
    float or ndarray
        Glide ratio [-].

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> float(glide_ratio(0.8, 0.05))
    16.0
    """
    _finite(Cl=Cl)
    _positive(Cd=Cd)
    return np.asarray(Cl) / Cd


def best_glide_cl(Cd0, AR, e=1.0):
    r"""Lift coefficient of maximum :math:`L/D` for a parabolic polar, :math:`\sqrt{\pi e AR C_{D0}}`.

    At this point induced drag equals parasite drag, so
    :math:`(L/D)_{\max} = \tfrac12\sqrt{\pi e AR / C_{D0}}`.

    Parameters
    ----------
    Cd0 : float or array_like
        Zero-lift drag coefficient [-], positive.
    AR : float or array_like
        Aspect ratio [-], positive.
    e : float or array_like, optional
        Oswald efficiency factor [-] in ``(0, 1]``.

    Returns
    -------
    float or ndarray
        Lift coefficient for best glide [-].

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(best_glide_cl(0.02, 8.0, 0.8)), 4)
    0.6341
    """
    _positive(Cd0=Cd0, AR=AR)
    _span_efficiency(e)
    return np.sqrt(np.pi * e * AR * np.asarray(Cd0))


def breguet_range(V, tsfc, lift_to_drag, W_initial, W_final, g=STANDARD_GRAVITY):
    r"""Jet-aircraft cruise range :math:`R = \dfrac{V}{g\,c_T}\dfrac{L}{D}\ln\dfrac{W_i}{W_f}` (Breguet).

    Integrates :math:`dR = V\,dt` while fuel burns at :math:`\dot m = c_T T`
    and thrust equals drag, :math:`T = W/(L/D)`, at constant speed and
    lift-to-drag ratio.

    Parameters
    ----------
    V : float or array_like
        Cruise true airspeed [m/s], positive.
    tsfc : float or array_like
        Thrust-specific fuel consumption [kg/(N s)], positive
        (e.g. ``1.6e-5`` for a modern turbofan).
    lift_to_drag : float or array_like
        Cruise lift-to-drag ratio [-], positive.
    W_initial, W_final : float or array_like
        Weight at the start and end of cruise [N], with ``W_initial >= W_final > 0``.
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    float or ndarray
        Still-air range [m].

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Aircraft Performance and Design* (McGraw-Hill, 1999),
    sec. 5.13, for the derivation of the range equation attributed to
    Louis Breguet.

    Examples
    --------
    >>> round(float(breguet_range(230.0, 1.6e-5, 18.0, 7.0e5, 5.0e5)) / 1e3)
    8878
    """
    _positive(V=V, tsfc=tsfc, lift_to_drag=lift_to_drag, W_final=W_final, g=g)
    _finite(W_initial=W_initial)
    if np.any(np.asarray(W_initial) < W_final):
        raise ValueError("W_initial must be at least W_final")
    return np.asarray(V) / (g * tsfc) * lift_to_drag * np.log(np.asarray(W_initial) / W_final)
