r"""Couette flow: viscous flow between moving walls and the rotating-cylinder viscometer.

Fluid between two parallel plates, one sliding at speed :math:`U`, takes a
linear velocity profile, and the shear stress :math:`\mu U/h` is the same
everywhere. A pressure gradient adds Poiseuille's parabola to the line.
Couette (1890) measured the viscosity of water this way, in the gap between
a fixed inner cylinder and a rotating outer one, from the torque on the
inner cylinder.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = ["couette_velocity", "couette_flow_rate", "couette_shear_stress", "couette_viscometer_torque"]


def couette_velocity(y, h, U, mu, dp_dx=0.0):
    r"""Velocity profile of plane Couette-Poiseuille flow.

    .. math::

        u(y) = U\frac{y}{h} - \frac{1}{2\mu}\frac{dp}{dx}\,y\,(h - y)

    Parameters
    ----------
    y : float or array_like
        Distance from the fixed plate [m], in [0, ``h``].
    h : float
        Gap between the plates [m], positive.
    U : float
        Speed of the moving plate at ``y = h`` [m/s].
    mu : float
        Dynamic viscosity [Pa s], positive.
    dp_dx : float, optional
        Pressure gradient along the flow [Pa/m]; negative pushes fluid
        forward. Zero gives pure Couette flow.

    Returns
    -------
    float or ndarray
        Velocity [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Steady, fully developed laminar flow
    of a Newtonian fluid: an exact solution of the Navier-Stokes
    equations. Backflow appears at the fixed wall when
    :math:`dp/dx > 2\mu U/h^2`.

    Examples
    --------
    >>> float(couette_velocity(0.5e-3, 1e-3, 2.0, 0.1))
    1.0
    """
    _finite(y=y, h=h, U=U, mu=mu, dp_dx=dp_dx)
    _positive(h=h, mu=mu)
    _interval("y", y, 0, h)
    y = np.asarray(y, dtype=float)
    return U * y / h - dp_dx / (2 * mu) * y * (h - y)


def couette_flow_rate(h, U, mu, dp_dx=0.0):
    r"""Flow per unit width between the plates, :math:`q = Uh/2 - h^3\,(dp/dx)/(12\mu)`.

    Parameters
    ----------
    h : float or array_like
        Gap [m], positive.
    U : float or array_like
        Speed of the moving plate [m/s].
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.
    dp_dx : float or array_like, optional
        Pressure gradient [Pa/m].

    Returns
    -------
    float or ndarray
        Volume flow per unit width [m^2/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The integral of
    :func:`couette_velocity` across the gap; the basis of Reynolds's
    lubrication theory of slider bearings.

    Examples
    --------
    >>> couette_flow_rate(1e-3, 2.0, 0.1)
    0.001
    """
    _finite(h=h, U=U, mu=mu, dp_dx=dp_dx)
    _positive(h=h, mu=mu)
    return U * h / 2 - h**3 * dp_dx / (12 * mu)


def couette_shear_stress(mu, U, h):
    r"""Shear stress in pure Couette flow, :math:`\tau = \mu U/h`.

    Parameters
    ----------
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.
    U : float or array_like
        Speed of the moving plate [m/s].
    h : float or array_like
        Gap [m], positive.

    Returns
    -------
    float or ndarray
        Shear stress [Pa], the same across the gap and on both plates.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Newton's law of viscosity applied to
    the linear profile, so a measured force per area and speed give
    :math:`\mu`.

    Examples
    --------
    >>> couette_shear_stress(0.1, 2.0, 1e-3)
    200.0
    """
    _finite(mu=mu, U=U, h=h)
    _positive(mu=mu, h=h)
    return mu * U / h


def couette_viscometer_torque(mu, omega, R_inner, R_outer, length):
    r"""Torque on the cylinders of a Couette viscometer.

    .. math::

        T = \frac{4\pi\mu\Omega L\,R_1^2R_2^2}{R_2^2 - R_1^2}

    Parameters
    ----------
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.
    omega : float or array_like
        Angular speed of one cylinder relative to the other [rad/s].
    R_inner : float or array_like
        Inner cylinder radius :math:`R_1` [m], positive.
    R_outer : float or array_like
        Outer cylinder radius :math:`R_2` [m], greater than ``R_inner``.
    length : float or array_like
        Wetted length of the cylinders [m], positive.

    Returns
    -------
    float or ndarray
        Torque [N m], the same on both cylinders.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The exact circular Couette solution
    :math:`v_\theta = Ar + B/r` with end effects neglected. For a narrow
    gap it tends to the plane result, shear stress :math:`\mu\Omega R/h`
    on area :math:`2\pi RL` at arm :math:`R`. Taylor vortices appear when
    the inner cylinder rotates fast; Couette rotated the outer one, which
    keeps the flow stable.

    References
    ----------
    M. Couette, "Études sur le frottement des liquides," *Ann. Chim.
    Phys.*, 6th ser., 21 (1890) 433-510.

    Examples
    --------
    >>> round(float(couette_viscometer_torque(1.0, 1.0, 1.0, np.sqrt(2.0), 1.0 / np.pi)), 6)
    8.0
    """
    _finite(mu=mu, omega=omega, R_inner=R_inner, R_outer=R_outer, length=length)
    _positive(mu=mu, R_inner=R_inner, R_outer=R_outer, length=length)
    if np.any(np.asarray(R_outer) <= np.asarray(R_inner)):
        raise ValueError("R_outer must be greater than R_inner")
    r1, r2 = R_inner * R_inner, R_outer * R_outer
    return 4 * np.pi * mu * omega * length * r1 * r2 / (r2 - r1)
