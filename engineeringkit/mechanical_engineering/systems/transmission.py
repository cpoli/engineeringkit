r"""Power transmission: gear trains and belt drives.

Euler showed that involute teeth keep the speed ratio of two gears constant,
so a gear pair is fixed by its tooth counts. Euler, and later Eytelwein for
engineers, found that friction lets a belt's tension grow exponentially with
the angle of wrap, :math:`T_1/T_2 = e^{\mu\theta}`. Lewis treated a gear
tooth as a cantilever of uniform strength inscribed in its profile, which
gives the bending stress :math:`\sigma = W_t/(F m Y)`, and Barth's factor
raises it for the dynamic load of a gear running at speed.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "gear_ratio",
    "gear_output_speed",
    "gear_output_torque",
    "belt_speed",
    "belt_tension_ratio",
    "lewis_form_factor",
    "lewis_bending_stress",
    "barth_velocity_factor",
]

# Lewis form factor of 20-degree full-depth spur teeth, Shigley's Mechanical
# Engineering Design, Table 14-2.
_LEWIS_TEETH = np.array([12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 26, 28, 30, 34, 38, 43, 50, 60, 75, 100, 150, 300, 400])
_LEWIS_Y = np.array(
    [
        0.245, 0.261, 0.277, 0.290, 0.296, 0.303, 0.309, 0.314, 0.322, 0.328, 0.331, 0.337, 0.346,
        0.353, 0.359, 0.371, 0.384, 0.397, 0.409, 0.422, 0.435, 0.447, 0.460, 0.472, 0.480,
    ]
)  # fmt: skip
_BARTH = {
    "cast": lambda V: (3.05 + V) / 3.05,
    "cut": lambda V: (6.1 + V) / 6.1,
    "hobbed": lambda V: (3.56 + np.sqrt(V)) / 3.56,
    "ground": lambda V: np.sqrt((5.56 + np.sqrt(V)) / 5.56),
}


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


def lewis_form_factor(teeth):
    r"""Lewis form factor :math:`Y` of a 20-degree full-depth involute spur tooth.

    Parameters
    ----------
    teeth : float or array_like
        Number of teeth [-], from 12 to 400; values between the tabulated
        tooth counts are interpolated linearly.

    Returns
    -------
    float or ndarray
        :math:`Y` [-], from 0.245 at 12 teeth to 0.480 at 400 (0.485 for a
        rack).

    Notes
    -----
    :math:`Y = 2xP/3` in Lewis's terms, with :math:`x` the depth from the
    root to where the inscribed parabola of uniform strength meets the tooth
    centreline and :math:`P` the diametral pitch (here :math:`1/m`).
    Values are those of Shigley's Table 14-2. Inputs are raw SI numbers with
    no dimensional metadata; see :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> float(lewis_form_factor(20))
    0.322
    """
    _interval("teeth", teeth, 12, 400)
    return np.interp(teeth, _LEWIS_TEETH, _LEWIS_Y)


def lewis_bending_stress(tangential_load, face_width, module, form_factor, velocity_factor=1.0):
    r"""Lewis bending stress at the root of a spur-gear tooth, :math:`\sigma = K_v W_t/(F m Y)`.

    Parameters
    ----------
    tangential_load : float or array_like
        Transmitted load :math:`W_t` at the pitch circle [N], nonnegative.
    face_width : float or array_like
        Face width :math:`F` [m], positive.
    module : float or array_like
        Module :math:`m` = pitch diameter over teeth [m], positive.
    form_factor : float or array_like
        Lewis form factor :math:`Y` [-], positive (see
        :func:`lewis_form_factor`).
    velocity_factor : float or array_like, optional
        Dynamic factor :math:`K_v` [-], at least one (see
        :func:`barth_velocity_factor`).

    Returns
    -------
    float or ndarray
        Bending stress [Pa].

    Notes
    -----
    The whole load is assumed to act at the tip of one tooth, with no stress
    concentration at the fillet and no radial compression; the AGMA method
    refines all three. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/mechanical_engineering`.

    References
    ----------
    W. Lewis, "Investigation of the strength of gear teeth," *Proc.
    Engineers' Club of Philadelphia* (1892) 16-23.

    Examples
    --------
    >>> round(float(lewis_bending_stress(1000.0, 0.03, 0.003, 0.322)) / 1e6, 3)
    34.507
    """
    _nonnegative(tangential_load=tangential_load)
    _positive(face_width=face_width, module=module, form_factor=form_factor)
    _interval("velocity_factor", velocity_factor, 1, np.inf, closed="left")
    return velocity_factor * np.asarray(tangential_load, dtype=float) / (face_width * module * form_factor)


def barth_velocity_factor(pitch_line_velocity, profile="cut"):
    r"""Barth's dynamic factor :math:`K_v` for spur gears running at speed.

    =========  =====================================
    profile    :math:`K_v`
    =========  =====================================
    cast       :math:`(3.05 + V)/3.05`
    cut        :math:`(6.1 + V)/6.1`
    hobbed     :math:`(3.56 + \sqrt V)/3.56`
    ground     :math:`\sqrt{(5.56 + \sqrt V)/5.56}`
    =========  =====================================

    Parameters
    ----------
    pitch_line_velocity : float or array_like
        :math:`V` [m/s], nonnegative.
    profile : {"cast", "cut", "hobbed", "ground"}, optional
        How the tooth profile was made.

    Returns
    -------
    float or ndarray
        :math:`K_v` [-], one at rest.

    Notes
    -----
    The constants are the SI forms of Barth's 1900-era empirical fits, as in
    Shigley's equation 14-4. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> float(barth_velocity_factor(6.1))
    2.0
    """
    _nonnegative(pitch_line_velocity=pitch_line_velocity)
    if profile not in _BARTH:
        raise ValueError(f"profile must be one of {sorted(_BARTH)}")
    return _BARTH[profile](np.asarray(pitch_line_velocity, dtype=float))
