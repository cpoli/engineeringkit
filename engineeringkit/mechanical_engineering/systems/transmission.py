r"""Power transmission: gear trains and belt drives.

Euler showed that involute teeth keep the speed ratio of two gears constant,
so a gear pair is fixed by its tooth counts. Euler, and later Eytelwein for
engineers, found that friction lets a belt's tension grow exponentially with
the angle of wrap, :math:`T_1/T_2 = e^{\mu\theta}`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["gear_ratio", "gear_output_speed", "gear_output_torque", "belt_speed", "belt_tension_ratio"]


def _teeth(teeth_driver, teeth_driven):
    _integer("teeth_driven", teeth_driven, 1)
    _integer("teeth_driver", teeth_driver, 1)


def gear_ratio(teeth_driven, teeth_driver):
    r"""Speed-reduction ratio :math:`N_\text{driven}/N_\text{driver}` of a gear pair.

    Parameters
    ----------
    teeth_driven : int
        Teeth on the driven gear, at least one.
    teeth_driver : int
        Teeth on the driving gear, at least one.

    Returns
    -------
    float
        Input speed over output speed [-]; above one for a reduction.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> gear_ratio(60, 20)
    3.0
    """
    _finite(teeth_driven=teeth_driven, teeth_driver=teeth_driver)
    _teeth(teeth_driver, teeth_driven)
    return teeth_driven / teeth_driver


def gear_output_speed(input_speed, teeth_driver, teeth_driven):
    r"""Output speed :math:`\omega_o = \omega_i N_i / N_o` of a gear pair.

    Parameters
    ----------
    input_speed : float or array_like
        Driver speed, in any angular-speed unit (rad/s or rpm).
    teeth_driver, teeth_driven : int
        Tooth counts, at least one.

    Returns
    -------
    float or ndarray
        Driven speed in the unit of ``input_speed``. The sign is preserved;
        the reversal of rotation in an external mesh is not applied.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> gear_output_speed(1500.0, 20, 60)
    500.0
    """
    _finite(input_speed=input_speed, teeth_driver=teeth_driver, teeth_driven=teeth_driven)
    _teeth(teeth_driver, teeth_driven)
    return input_speed * teeth_driver / teeth_driven


def gear_output_torque(input_torque, teeth_driver, teeth_driven, efficiency=1):
    r"""Output torque :math:`T_o = \eta T_i N_o / N_i` of a gear pair.

    Parameters
    ----------
    input_torque : float or array_like
        Driver torque [N m].
    teeth_driver, teeth_driven : int
        Tooth counts, at least one.
    efficiency : float, optional
        Mesh power efficiency [-] in [0, 1].

    Returns
    -------
    float or ndarray
        Driven torque [N m]. Output power is ``efficiency`` times input power.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> gear_output_torque(10.0, 20, 60, efficiency=0.9)
    27.0
    """
    _finite(input_torque=input_torque, teeth_driver=teeth_driver, teeth_driven=teeth_driven, efficiency=efficiency)
    _teeth(teeth_driver, teeth_driven)
    _interval("efficiency", efficiency, 0, 1)
    return input_torque * teeth_driven / teeth_driver * efficiency


def belt_speed(diameter, rpm):
    r"""Linear speed :math:`v = \pi D n / 60` of a belt on a pulley.

    Parameters
    ----------
    diameter : float or array_like
        Pitch diameter of the pulley [m], positive.
    rpm : float or array_like
        Pulley speed [rev/min].

    Returns
    -------
    float or ndarray
        Belt speed [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata, except the speed
    in rev/min; see :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(belt_speed(1.0, 60.0)), 6)
    3.141593
    """
    _finite(diameter=diameter, rpm=rpm)
    _positive(diameter=diameter)
    return np.pi * diameter * rpm / 60


def belt_tension_ratio(mu, wrap_angle, groove_angle=None):
    r"""Euler-Eytelwein ratio of tight-side to slack-side tension at the point of slip.

    .. math:: \frac{T_1}{T_2} = \exp\!\left(\frac{\mu\theta}{\sin(\beta/2)}\right)

    Parameters
    ----------
    mu : float or array_like
        Coulomb friction coefficient between belt and pulley [-],
        nonnegative.
    wrap_angle : float or array_like
        Angle of contact :math:`\theta` [rad], nonnegative.
    groove_angle : float, optional
        Included angle :math:`\beta` of a V-belt groove [rad] in (0, pi).
        Omit for a flat belt or rope on a drum (:math:`\sin(\beta/2) = 1`).

    Returns
    -------
    float or ndarray
        :math:`T_1/T_2` [-], at least one.

    Notes
    -----
    The belt is massless (no centrifugal tension) and inextensible. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    L. Euler, "Remarques sur l'effect du frottement dans l'équilibre,"
    *Mém. Acad. Sci. Berlin* 18 (1762, published 1769) 265-278;
    J. A. Eytelwein, *Handbuch der Statik fester Körper* (Berlin, 1808).

    Examples
    --------
    >>> round(float(belt_tension_ratio(0.3, np.pi)), 4)
    2.5663
    """
    _nonnegative(mu=mu, wrap_angle=wrap_angle)
    wedge = 1.0
    if groove_angle is not None:
        _interval("groove_angle", groove_angle, 0, np.pi, closed=False)
        wedge = np.sin(groove_angle / 2)
    return np.exp(np.asarray(mu) * wrap_angle / wedge)
