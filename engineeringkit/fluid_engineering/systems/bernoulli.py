r"""Bernoulli's equation and continuity for steady incompressible flow.

Along a streamline of steady, inviscid, incompressible flow the total head
:math:`p/(\rho g) + v^2/(2g) + z` is constant (Daniel Bernoulli 1738, Euler
1757), and continuity fixes the velocity in each cross-section,
:math:`A_1 v_1 = A_2 v_2`. Flow meters and the energy and hydraulic grade
lines of pipe systems are built on these two relations.
"""

from __future__ import annotations

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G

__all__ = ["bernoulli_pressure", "continuity_velocity", "dynamic_pressure", "velocity_head"]


def bernoulli_pressure(p1, rho, v1, v2, z1=0, z2=0, g=G):
    r"""Pressure at point 2 from Bernoulli's equation along a streamline.

    .. math::

        p_2 = p_1 + \tfrac12\rho(v_1^2 - v_2^2) + \rho g (z_1 - z_2)

    Parameters
    ----------
    p1 : float or array_like
        Pressure at point 1 [Pa].
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    v1, v2 : float or array_like
        Speeds at points 1 and 2 [m/s].
    z1, z2 : float or array_like, optional
        Elevations of points 1 and 2 [m].
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Pressure at point 2 [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Steady, incompressible, inviscid flow
    with no pump or turbine between the points; add losses separately.

    References
    ----------
    D. Bernoulli, *Hydrodynamica* (Strasbourg, 1738).

    Examples
    --------
    >>> bernoulli_pressure(2000.0, 1000.0, 0.0, 2.0)
    0.0
    """
    _finite(p1=p1, rho=rho, v1=v1, v2=v2, z1=z1, z2=z2, g=g)
    _positive(rho=rho, g=g)
    return p1 + 0.5 * rho * (v1 * v1 - v2 * v2) + rho * g * (z1 - z2)


def continuity_velocity(A1, v1, A2):
    r"""Mean velocity in section 2 from continuity, :math:`v_2 = A_1 v_1/A_2`.

    Parameters
    ----------
    A1 : float or array_like
        Area of section 1 [m^2], positive.
    v1 : float or array_like
        Mean velocity in section 1 [m/s].
    A2 : float or array_like
        Area of section 2 [m^2], positive.

    Returns
    -------
    float or ndarray
        Mean velocity in section 2 [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> continuity_velocity(4.0, 1.0, 1.0)
    4.0
    """
    _finite(A1=A1, v1=v1, A2=A2)
    _positive(A1=A1, A2=A2)
    return A1 * v1 / A2


def dynamic_pressure(rho, v):
    r"""Dynamic pressure :math:`q = \tfrac12\rho v^2`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    v : float or array_like
        Flow speed [m/s].

    Returns
    -------
    float or ndarray
        Dynamic pressure [Pa], the rise from static to stagnation pressure.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> dynamic_pressure(1000.0, 2.0)
    2000.0
    """
    _finite(rho=rho, v=v)
    _positive(rho=rho)
    return 0.5 * rho * v * v


def velocity_head(v, g=G):
    r"""Velocity head :math:`v^2/(2g)`.

    Parameters
    ----------
    v : float or array_like
        Mean flow speed [m/s].
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Velocity head [m], the gap between the energy and hydraulic grade lines.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> velocity_head(2.0, g=10.0)
    0.2
    """
    _finite(v=v, g=g)
    _positive(g=g)
    return v * v / (2 * g)
