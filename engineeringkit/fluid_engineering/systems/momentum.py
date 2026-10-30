r"""The momentum equation: forces of jets on vanes and of flow on pipe bends.

For steady flow through a control volume the net force on the fluid is the
momentum leaving per second minus the momentum entering,
:math:`\sum\mathbf F = \sum \dot m\,\mathbf V_{out} - \sum \dot m\,\mathbf V_{in}`.
A jet turned by a vane pushes the vane; flow turned by a bend, with the
pressure at each end, pushes the bend off its supports. Pelton's wheel
(1880) turns the jet almost back on itself to take nearly all its energy.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["jet_force", "bend_force", "pelton_wheel_power"]


def jet_force(rho, area, v, angle, vane_speed=0.0):
    r"""Force of a free jet on a single vane that turns it through ``angle``.

    .. math::

        F_x = \rho A (V - u)^2 (1 - \cos\theta),\qquad
        F_y = \rho A (V - u)^2 \sin\theta

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    area : float or array_like
        Jet cross-section [m^2], positive.
    v : float or array_like
        Jet speed [m/s], nonnegative.
    angle : float or array_like
        Angle through which the vane turns the jet [rad], in [0, π]:
        π/2 for a flat plate normal to the jet, π for a full reversal.
    vane_speed : float or array_like, optional
        Speed of the vane along the jet [m/s], in [0, ``v``].

    Returns
    -------
    tuple of float or ndarray
        ``(Fx, Fy)`` [N]: the force on the vane along the jet and
        across it, toward the side the jet is turned away from.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The jet is at atmospheric pressure
    throughout, gravity and vane friction are neglected, so the relative
    speed :math:`V - u` is unchanged over the vane. A single moving vane
    receives only :math:`\rho A (V - u)` per second because it runs away
    from the jet; a wheel of vanes receives the whole flow, see
    :func:`pelton_wheel_power`.

    Examples
    --------
    >>> round(float(jet_force(1000.0, 0.01, 10.0, np.pi / 2)[0]), 6)
    1000.0
    """
    _finite(rho=rho, area=area, v=v, angle=angle, vane_speed=vane_speed)
    _positive(rho=rho, area=area)
    _nonnegative(v=v, vane_speed=vane_speed)
    _interval("angle", angle, 0, np.pi)
    if np.any(np.asarray(vane_speed) > np.asarray(v)):
        raise ValueError("vane_speed must not exceed the jet speed")
    momentum = rho * area * (v - vane_speed) ** 2
    return momentum * (1 - np.cos(angle)), momentum * np.sin(angle)


def bend_force(rho, Q, A1, A2, p1, p2, angle):
    r"""Force of the flow on a horizontal reducing bend, which its anchor must resist.

    Flow enters along :math:`+x` and leaves at ``angle`` from it. The
    momentum balance on the fluid inside the bend gives the force on the
    bend:

    .. math::

        F_x = p_1A_1 + \rho QV_1 - (p_2A_2 + \rho QV_2)\cos\theta,\qquad
        F_y = -(p_2A_2 + \rho QV_2)\sin\theta

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    Q : float or array_like
        Volume flow [m^3/s], nonnegative.
    A1 : float or array_like
        Inlet area [m^2], positive.
    A2 : float or array_like
        Outlet area [m^2], positive.
    p1 : float or array_like
        Gauge pressure at the inlet [Pa].
    p2 : float or array_like
        Gauge pressure at the outlet [Pa]; from
        :func:`bernoulli_pressure` less any loss.
    angle : float or array_like
        Turning angle :math:`\theta` [rad], in [0, π], counterclockwise.

    Returns
    -------
    tuple of float or ndarray
        ``(Fx, Fy)`` [N]: the force of the fluid on the bend. The anchor
        supplies the opposite force.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Gauge pressures, because the
    atmosphere acts on the outside of the bend too. Uniform velocity at
    each end; the weight of the bend and its water acts vertically and is
    not included.

    Examples
    --------
    A 180° return bend with no pressure: twice the momentum flux.

    >>> round(float(bend_force(1000.0, 0.01, 0.01, 0.01, 0.0, 0.0, np.pi)[0]), 6)
    20.0
    """
    _finite(rho=rho, Q=Q, A1=A1, A2=A2, p1=p1, p2=p2, angle=angle)
    _positive(rho=rho, A1=A1, A2=A2)
    _nonnegative(Q=Q)
    _interval("angle", angle, 0, np.pi)
    inlet = p1 * A1 + rho * Q * Q / A1
    outlet = p2 * A2 + rho * Q * Q / A2
    return inlet - outlet * np.cos(angle), -outlet * np.sin(angle)


def pelton_wheel_power(rho, Q, v, bucket_speed, angle=np.radians(165)):
    r"""Power a Pelton wheel takes from its jet, :math:`P = \rho Q u (V - u)(1 - \cos\theta)`.

    Parameters
    ----------
    rho : float or array_like
        Water density [kg/m^3], positive.
    Q : float or array_like
        Jet volume flow [m^3/s], nonnegative.
    v : float or array_like
        Jet speed [m/s], nonnegative; :math:`C_v\sqrt{2gH}` from the head.
    bucket_speed : float or array_like
        Bucket speed at the pitch circle :math:`u = \omega r` [m/s], in
        [0, ``v``].
    angle : float or array_like, optional
        Angle through which the buckets turn the jet [rad], in [0, π];
        165° by default, short of 180° so the water clears the next bucket.

    Returns
    -------
    float or ndarray
        Power on the wheel [W]. It is greatest at :math:`u = V/2`, where
        the hydraulic efficiency is :math:`(1 - \cos\theta)/2`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Unlike one vane in
    :func:`jet_force`, a wheel always has a bucket in the jet, so the whole
    flow :math:`Q` is turned. Bucket friction is neglected.

    References
    ----------
    L. A. Pelton, "Water-wheel," US Patent 233,692 (26 October 1880).

    Examples
    --------
    >>> float(pelton_wheel_power(1000.0, 0.1, 20.0, 10.0, np.pi))
    20000.0
    """
    _finite(rho=rho, Q=Q, v=v, bucket_speed=bucket_speed, angle=angle)
    _positive(rho=rho)
    _nonnegative(Q=Q, v=v, bucket_speed=bucket_speed)
    _interval("angle", angle, 0, np.pi)
    if np.any(np.asarray(bucket_speed) > np.asarray(v)):
        raise ValueError("bucket_speed must not exceed the jet speed")
    return rho * Q * bucket_speed * (v - bucket_speed) * (1 - np.cos(angle))
