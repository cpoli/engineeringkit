r"""Rigid-body dynamics: inertia, work, energy, power and momentum.

Huygens found the centre of oscillation of a compound pendulum from its
moment of inertia, and Steiner stated the parallel-axis theorem
:math:`I = I_G + md^2` in general. Coriolis gave *work* and *kinetic
energy*, with the factor one half, the names and meanings engineers still
use.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "kinetic_energy",
    "rotational_energy",
    "work_constant_force",
    "power_torque",
    "angular_momentum",
    "centripetal_force",
    "flywheel_energy",
    "disk_inertia",
    "solid_sphere_inertia",
    "parallel_axis",
]


def kinetic_energy(m, v):
    r"""Translational kinetic energy :math:`T = \tfrac12 m v^2`.

    Parameters
    ----------
    m : float or array_like
        Mass [kg], positive.
    v : float or array_like
        Speed [m/s].

    Returns
    -------
    float or ndarray
        Kinetic energy [J].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    G. Coriolis, *Du calcul de l'effet des machines* (Carilian-Goeury,
    Paris, 1829).

    Examples
    --------
    >>> kinetic_energy(2.0, 3.0)
    9.0
    """
    _finite(m=m, v=v)
    _positive(m=m)
    return 0.5 * m * v * v


def rotational_energy(I, omega):
    r"""Kinetic energy of rotation about a fixed axis, :math:`\tfrac12 I\omega^2`.

    Parameters
    ----------
    I : float or array_like
        Moment of inertia about the axis [kg m^2], positive.
    omega : float or array_like
        Angular speed [rad/s].

    Returns
    -------
    float or ndarray
        Kinetic energy [J].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> rotational_energy(2.0, 3.0)
    9.0
    """
    _finite(I=I, omega=omega)
    _positive(I=I)
    return 0.5 * I * omega * omega


def work_constant_force(F, d, theta=0):
    r"""Work :math:`W = F d \cos\theta` of a constant force over a straight displacement.

    Parameters
    ----------
    F : float or array_like
        Force magnitude [N].
    d : float or array_like
        Displacement magnitude [m].
    theta : float or array_like, optional
        Angle between force and displacement [rad].

    Returns
    -------
    float or ndarray
        Work done by the force [J]; negative when it opposes the motion.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(work_constant_force(10.0, 2.0, np.pi / 3)), 12)
    10.0
    """
    _finite(F=F, d=d, theta=theta)
    return F * d * np.cos(theta)


def power_torque(torque, omega):
    r"""Shaft power :math:`P = T\omega`.

    Parameters
    ----------
    torque : float or array_like
        Torque [N m].
    omega : float or array_like
        Angular speed [rad/s].

    Returns
    -------
    float or ndarray
        Power [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> power_torque(10.0, 100.0)
    1000.0
    """
    _finite(torque=torque, omega=omega)
    return torque * omega


def angular_momentum(I, omega):
    r"""Angular momentum :math:`H = I\omega` about a fixed axis.

    Parameters
    ----------
    I : float or array_like
        Moment of inertia [kg m^2], positive.
    omega : float or array_like
        Angular speed [rad/s].

    Returns
    -------
    float or ndarray
        Angular momentum [kg m^2/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> angular_momentum(2.0, 3.0)
    6.0
    """
    _finite(I=I, omega=omega)
    _positive(I=I)
    return I * omega


def centripetal_force(m, v, r):
    r"""Force :math:`m v^2 / r` that holds a mass on a circular path.

    Parameters
    ----------
    m : float or array_like
        Mass [kg], positive.
    v : float or array_like
        Tangential speed [m/s].
    r : float or array_like
        Path radius [m], positive.

    Returns
    -------
    float or ndarray
        Force directed toward the centre [N].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> centripetal_force(2.0, 3.0, 1.0)
    18.0
    """
    _finite(m=m, v=v, r=r)
    _positive(m=m, r=r)
    return m * v * v / r


def flywheel_energy(I, w1, w2):
    r"""Energy a flywheel releases slowing from :math:`\omega_1` to :math:`\omega_2`.

    .. math:: \Delta E = \tfrac12 I (\omega_1^2 - \omega_2^2)

    Parameters
    ----------
    I : float or array_like
        Moment of inertia [kg m^2], positive.
    w1, w2 : float or array_like
        Initial and final angular speeds [rad/s].

    Returns
    -------
    float or ndarray
        Energy released [J]; negative if the flywheel speeds up.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> flywheel_energy(2.0, 3.0, 1.0)
    8.0
    """
    _finite(I=I, w1=w1, w2=w2)
    _positive(I=I)
    return 0.5 * I * (w1 * w1 - w2 * w2)


def disk_inertia(m, r):
    r"""Moment of inertia :math:`\tfrac12 m r^2` of a uniform disk about its axis.

    Parameters
    ----------
    m : float or array_like
        Mass [kg], positive.
    r : float or array_like
        Radius [m], positive.

    Returns
    -------
    float or ndarray
        Moment of inertia [kg m^2].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> disk_inertia(2.0, 1.0)
    1.0
    """
    _finite(m=m, r=r)
    _positive(m=m, r=r)
    return 0.5 * m * r * r


def solid_sphere_inertia(m, r):
    r"""Moment of inertia :math:`\tfrac25 m r^2` of a uniform sphere about a diameter.

    Parameters
    ----------
    m : float or array_like
        Mass [kg], positive.
    r : float or array_like
        Radius [m], positive.

    Returns
    -------
    float or ndarray
        Moment of inertia [kg m^2].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> solid_sphere_inertia(5.0, 1.0)
    2.0
    """
    _finite(m=m, r=r)
    _positive(m=m, r=r)
    return 0.4 * m * r * r


def parallel_axis(I_centroid, m, d):
    r"""Steiner's parallel-axis theorem, :math:`I = I_G + m d^2`.

    Parameters
    ----------
    I_centroid : float or array_like
        Moment of inertia about an axis through the centre of mass
        [kg m^2], nonnegative.
    m : float or array_like
        Mass [kg], positive.
    d : float or array_like
        Distance from that axis to the parallel axis [m].

    Returns
    -------
    float or ndarray
        Moment of inertia about the parallel axis [kg m^2].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    C. Huygens, *Horologium oscillatorium* (Muguet, Paris, 1673), part IV;
    J. Steiner, "Von dem Krümmungs-Schwerpuncte ebener Curven," *J. Reine
    Angew. Math.* 21 (1840) 33-63.

    Examples
    --------
    >>> parallel_axis(1.0, 2.0, 1.0)
    3.0
    """
    _finite(I_centroid=I_centroid, m=m, d=d)
    _nonnegative(I_centroid=I_centroid)
    _positive(m=m)
    return I_centroid + m * d * d
