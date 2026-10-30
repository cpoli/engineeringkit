r"""Drag on a sphere: Stokes's creeping-flow law and the drag curve.

Stokes (1851) solved the slow flow past a sphere and found the drag
:math:`3\pi\mu d V`, so :math:`C_D = 24/Re`. It holds while inertia is
negligible, :math:`Re \lesssim 1`. At higher Reynolds numbers a wake forms
and :math:`C_D` levels off near 0.4 until the boundary layer turns
turbulent at about :math:`Re = 2\times10^5` (the drag crisis). Setting the
drag equal to the weight less the buoyancy gives the terminal velocity of a
falling particle or a rising bubble.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G

__all__ = ["stokes_drag", "stokes_settling_velocity", "sphere_drag_coefficient", "terminal_velocity"]

DRAG_CRISIS_RE = 2e5
"""Reynolds number above which the sphere drag curve is not fitted [-]."""


def stokes_drag(mu, v, d):
    r"""Stokes's drag on a sphere in creeping flow, :math:`F = 3\pi\mu d V`.

    Parameters
    ----------
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.
    v : float or array_like
        Speed of the sphere relative to the fluid [m/s]; the force has its
        sign.
    d : float or array_like
        Sphere diameter [m], positive.

    Returns
    -------
    float or ndarray
        Drag force [N], opposing the motion.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Valid for :math:`Re = \rho V d/\mu`
    below about 1. Two thirds of the drag is skin friction and one third
    pressure.

    References
    ----------
    G. G. Stokes, "On the effect of the internal friction of fluids on the
    motion of pendulums," *Trans. Camb. Phil. Soc.* 9, Part II (1851)
    8-106 (read 9 December 1850).

    Examples
    --------
    >>> round(float(stokes_drag(1e-3, 0.01, 1e-3 / np.pi)), 12)
    3e-08
    """
    _finite(mu=mu, v=v, d=d)
    _positive(mu=mu, d=d)
    return 3 * np.pi * mu * d * v


def stokes_settling_velocity(d, rho_particle, rho, mu, g=G):
    r"""Terminal velocity of a small sphere in Stokes flow, :math:`V = (\rho_p - \rho) g d^2/(18\mu)`.

    Parameters
    ----------
    d : float or array_like
        Sphere diameter [m], positive.
    rho_particle : float or array_like
        Sphere density [kg/m^3], positive.
    rho : float or array_like
        Fluid density [kg/m^3], nonnegative.
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Settling velocity [m/s], downward; negative when the sphere is
        lighter than the fluid and rises.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Check that the Reynolds number of
    the result is below about 1; otherwise use :func:`terminal_velocity`.

    Examples
    --------
    >>> round(float(stokes_settling_velocity(1e-4, 2650.0, 1000.0, 1e-3, g=10.0)), 6)
    0.009167
    """
    _finite(d=d, rho_particle=rho_particle, rho=rho, mu=mu, g=g)
    _positive(d=d, rho_particle=rho_particle, mu=mu, g=g)
    _nonnegative(rho=rho)
    return (rho_particle - rho) * g * d * d / (18 * mu)


def sphere_drag_coefficient(Re):
    r"""Drag coefficient of a smooth sphere, :math:`C_D \approx 24/Re + 6/(1 + \sqrt{Re}) + 0.4`.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number :math:`\rho V d/\mu` [-], in (0, 2e5].

    Returns
    -------
    float or ndarray
        Drag coefficient on the frontal area :math:`\pi d^2/4` [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. White's curve fit, within about 10%
    of measurements up to the drag crisis; it tends to Stokes's
    :math:`24/Re` as :math:`Re \to 0`. Above :math:`2\times10^5` the
    boundary layer turns turbulent and :math:`C_D` drops to about 0.1,
    which this fit does not represent, so those values are rejected.

    References
    ----------
    F. M. White, *Viscous Fluid Flow*, 3rd ed. (New York: McGraw-Hill,
    2006), chapter 3.

    Examples
    --------
    >>> float(sphere_drag_coefficient(1.0))
    27.4
    """
    _finite(Re=Re)
    _positive(Re=Re)
    if np.any(np.asarray(Re) > DRAG_CRISIS_RE):
        raise ValueError(f"Re must not exceed {DRAG_CRISIS_RE:g}, the drag crisis")
    return 24 / Re + 6 / (1 + np.sqrt(Re)) + 0.4


def terminal_velocity(d, rho_particle, rho, mu, g=G):
    r"""Terminal velocity of a sphere falling through a fluid, from the drag curve.

    Drag balances weight less buoyancy when
    :math:`C_D(Re)\,Re^2 = \tfrac43\,(\rho_p - \rho)\rho g d^3/\mu^2`; the
    left side, with :func:`sphere_drag_coefficient`, grows steadily with
    :math:`Re`, so Brent's method finds the single root.

    Parameters
    ----------
    d : float
        Sphere diameter [m], positive.
    rho_particle : float
        Sphere density [kg/m^3], greater than ``rho``.
    rho : float
        Fluid density [kg/m^3], positive.
    mu : float
        Dynamic viscosity [Pa s], positive.
    g : float, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float
        Terminal velocity [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Scalar inputs only. Raises
    ``ValueError`` when the answer lies beyond the drag crisis
    (:math:`Re > 2\times10^5`). Equal to :func:`stokes_settling_velocity`
    in the limit of small spheres.

    Examples
    --------
    A 1 mm sand grain in water:

    >>> round(terminal_velocity(1e-3, 2650.0, 1000.0, 1e-3), 3)
    0.145
    """
    _finite(d=d, rho_particle=rho_particle, rho=rho, mu=mu, g=g)
    _positive(d=d, rho=rho, mu=mu, g=g)
    if rho_particle <= rho:
        raise ValueError("rho_particle must exceed rho for the sphere to fall")
    best = 4 / 3 * (rho_particle - rho) * rho * g * d**3 / mu**2

    def residual(Re):
        # C_D Re^2, written so that it is defined at Re = 0.
        return 24 * Re + 6 * Re * Re / (1 + np.sqrt(Re)) + 0.4 * Re * Re - best

    if residual(DRAG_CRISIS_RE) < 0:
        raise ValueError("terminal velocity lies beyond the drag crisis, Re > 2e5")
    Re = optimize.brentq(residual, 0.0, DRAG_CRISIS_RE, xtol=1e-12, rtol=1e-14)
    return float(Re * mu / (rho * d))
