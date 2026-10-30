r"""Bolted joints: preload, bolt and member stiffness and the load factor.

A tightened bolt is a stretched spring clamping compressed members. An
external load :math:`P` that tries to separate them stretches the bolt a
little more and relieves the members by the same amount, so the two share
:math:`P` in proportion to their stiffnesses: the bolt takes :math:`CP`,
with :math:`C = k_b/(k_b + k_m)`, until the members part at
:math:`P_0 = F_i/(1 - C)`. Rötscher drew this as the joint diagram. The
members' stiffness comes from the cones of compressed material spreading out
from under the bolt head and nut.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.mechanical_engineering.core.base import BoltedJoint

__all__ = ["bolt_stiffness", "member_stiffness", "bolted_joint", "preload_torque"]


def bolt_stiffness(E, shank_area, shank_length, thread_area, thread_length):
    r"""Axial stiffness of a bolt's unthreaded shank and threaded length in series.

    .. math:: k_b = \frac{A_d A_t E}{A_d l_t + A_t l_d}

    Parameters
    ----------
    E : float
        Young's modulus of the bolt [Pa], positive.
    shank_area : float
        Major-diameter area :math:`A_d` [m^2], positive.
    shank_length : float
        Unthreaded length within the grip :math:`l_d` [m], nonnegative.
    thread_area : float
        Tensile-stress area :math:`A_t` [m^2], positive.
    thread_length : float
        Threaded length within the grip :math:`l_t` [m], nonnegative.

    Returns
    -------
    float
        :math:`k_b` [N/m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(bolt_stiffness(200e9, 1e-4, 0.0, 1e-4, 0.05) / 1e6, 6)
    400.0
    """
    _positive(E=E, shank_area=shank_area, thread_area=thread_area)
    _nonnegative(shank_length=shank_length, thread_length=thread_length)
    if shank_length + thread_length == 0:
        raise ValueError("the grip length shank_length + thread_length must be positive")
    return float(shank_area * thread_area * E / (shank_area * thread_length + thread_area * shank_length))


def member_stiffness(E, d, grip, washer_diameter=None, cone_angle=np.pi / 6):
    r"""Stiffness of two identical clamped members from the pressure-cone (frustum) model.

    Each member is a hollow frustum of half-angle :math:`\alpha` spreading
    from the washer face :math:`D` round the bolt hole :math:`d`; the two,
    each of thickness :math:`l/2`, act in series:

    .. math::

        k_m = \frac{\pi E d \tan\alpha}
        {2\ln\dfrac{(l\tan\alpha + D - d)(D + d)}{(l\tan\alpha + D + d)(D - d)}}

    Parameters
    ----------
    E : float
        Young's modulus of the members [Pa], positive.
    d : float
        Bolt diameter [m], positive.
    grip : float
        Total clamped thickness :math:`l` [m], positive.
    washer_diameter : float, optional
        Washer-face diameter :math:`D` [m], greater than ``d``; ``1.5 d`` if
        omitted.
    cone_angle : float, optional
        Frustum half-angle :math:`\alpha` [rad] in (0, pi/2); 30 degrees is
        the usual choice.

    Returns
    -------
    float
        :math:`k_m` [N/m].

    Notes
    -----
    With :math:`D = 1.5d` and 30 degrees this is Shigley's equation 8-22.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    F. Rötscher, *Die Maschinenelemente* (Springer, Berlin, 1927).

    Examples
    --------
    >>> round(member_stiffness(200e9, 0.012, 0.04) / 1e9, 4)
    2.1595
    """
    _positive(E=E, d=d, grip=grip)
    D = 1.5 * d if washer_diameter is None else washer_diameter
    _positive(washer_diameter=D)
    if D <= d:
        raise ValueError("washer_diameter must exceed the bolt diameter d")
    _interval("cone_angle", cone_angle, 0, np.pi / 2, closed=False)
    t = grip * np.tan(cone_angle)
    return float(np.pi * E * d * np.tan(cone_angle) / (2 * np.log((t + D - d) * (D + d) / ((t + D + d) * (D - d)))))


def bolted_joint(bolt_stiffness, member_stiffness, preload, load) -> BoltedJoint:
    r"""Share an external separating load between a preloaded bolt and its members.

    Parameters
    ----------
    bolt_stiffness : float
        :math:`k_b` [N/m], positive (see :func:`bolt_stiffness`).
    member_stiffness : float
        :math:`k_m` [N/m], positive (see :func:`member_stiffness`).
    preload : float
        Initial tension :math:`F_i` [N], nonnegative.
    load : float
        External tensile load on the joint :math:`P` [N], positive.

    Returns
    -------
    BoltedJoint
        Joint constant, bolt and member loads, separation load and
        ``load_factor(proof_load)``.

    Notes
    -----
    The members are assumed to stay in contact, and the load to act at the
    interfaces under the head and nut, until separation. Inputs are raw SI
    numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> joint = bolted_joint(1e9, 3e9, 20e3, 8e3)
    >>> joint.joint_constant, joint.bolt_load, joint.member_load
    (0.25, 22000.0, 14000.0)
    """
    _positive(bolt_stiffness=bolt_stiffness, member_stiffness=member_stiffness, load=load)
    _nonnegative(preload=preload)
    return BoltedJoint(float(bolt_stiffness), float(member_stiffness), float(preload), float(load))


def preload_torque(preload, d, torque_coefficient=0.2):
    r"""Wrench torque to reach a bolt preload, :math:`T = K F_i d`.

    Parameters
    ----------
    preload : float or array_like
        :math:`F_i` [N], nonnegative.
    d : float or array_like
        Nominal bolt diameter [m], positive.
    torque_coefficient : float or array_like, optional
        :math:`K` [-], positive; about 0.2 for unlubricated steel and 0.15
        to 0.18 when lubricated.

    Returns
    -------
    float or ndarray
        Torque [N m].

    Notes
    -----
    :math:`K` lumps thread and collar friction and scatters by about
    :math:`\pm 25\%`, so torque control sets the preload only roughly.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(preload_torque(20e3, 0.012)), 6)
    48.0
    """
    _nonnegative(preload=preload)
    _positive(d=d, torque_coefficient=torque_coefficient)
    return torque_coefficient * np.asarray(preload, dtype=float) * d
