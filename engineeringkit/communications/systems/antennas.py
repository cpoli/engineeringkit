r"""Antennas: gain and effective aperture, directivity from a pattern, and the half-wave dipole.

Friis wrote his transmission formula in terms of the antennas' effective
areas, related to gain by :math:`A_e = G\lambda^2/4\pi`. Carter's
induced-EMF method (1932) gave the half-wave dipole its radiation
resistance of about 73 ohms and, with its pattern
:math:`\cos(\tfrac\pi2\cos\theta)/\sin\theta`, a directivity of 1.64
(2.15 dBi).
"""

from __future__ import annotations

import numpy as np
from scipy import special
from scipy.integrate import trapezoid

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.communications.systems.propagation import SPEED_OF_LIGHT

__all__ = [
    "FREE_SPACE_IMPEDANCE",
    "effective_aperture",
    "aperture_gain",
    "directivity",
    "half_wave_dipole_pattern",
    "half_wave_dipole_directivity",
    "half_wave_dipole_resistance",
]

FREE_SPACE_IMPEDANCE = 376.730313668
"""Impedance of free space :math:`\\eta_0 = \\mu_0 c` [ohm] (CODATA 2018)."""


def _cin_2pi():
    # Cin(x) = gamma + ln x - Ci(x), the integral of (1 - cos t)/t from 0 to x.
    return float(np.euler_gamma + np.log(2 * np.pi) - special.sici(2 * np.pi)[1])


def effective_aperture(gain, f, c=SPEED_OF_LIGHT):
    r"""Effective aperture of an antenna, :math:`A_e = G\lambda^2/4\pi`.

    Parameters
    ----------
    gain : float or array_like
        Antenna gain [-] (linear, not dBi), positive.
    f : float or array_like
        Frequency [Hz], positive.
    c : float, optional
        Propagation speed [m/s], positive.

    Returns
    -------
    float or ndarray
        Effective area [m^2]: the received power is the incident power
        density times this area.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    H. T. Friis, "A note on a simple transmission formula," *Proc. IRE* 34
    (1946) 254-256.

    Examples
    --------
    >>> round(float(effective_aperture(1.0, 299792458.0)), 6)
    0.079577
    """
    _positive(gain=gain, f=f, c=c)
    return np.asarray(gain, dtype=float) * (c / np.asarray(f, dtype=float)) ** 2 / (4 * np.pi)


def aperture_gain(area, f, efficiency=1.0, c=SPEED_OF_LIGHT):
    r"""Gain of an aperture antenna such as a dish, :math:`G = 4\pi\eta_a A/\lambda^2`.

    Parameters
    ----------
    area : float or array_like
        Physical aperture area [m^2], positive.
    f : float or array_like
        Frequency [Hz], positive.
    efficiency : float or array_like, optional
        Aperture efficiency :math:`\eta_a` [-] in ``(0, 1]``; 0.55 to 0.7 is
        typical of a parabolic dish.
    c : float, optional
        Propagation speed [m/s], positive.

    Returns
    -------
    float or ndarray
        Gain [-] (linear).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(10 * np.log10(aperture_gain(np.pi * 1.5**2, 12e9, 0.6))), 2)
    49.31
    """
    _positive(area=area, f=f, c=c)
    _interval("efficiency", efficiency, 0, 1, closed="right")
    wavelength = c / np.asarray(f, dtype=float)
    return 4 * np.pi * np.asarray(efficiency, dtype=float) * np.asarray(area, dtype=float) / wavelength**2


def directivity(theta, power_pattern):
    r"""Directivity of an azimuthally symmetric pattern, :math:`D = 2U_{max}\big/\int_0^\pi U(\theta)\sin\theta\,d\theta`.

    Parameters
    ----------
    theta : array_like
        Increasing polar angles [rad] covering ``[0, pi]``.
    power_pattern : array_like
        Radiation intensity :math:`U(\theta)` at those angles [any unit],
        nonnegative, not all zero.

    Returns
    -------
    float
        Directivity [-] (linear), by the trapezoidal rule; 1.5 for a short
        (Hertzian) dipole, whose pattern is :math:`\sin^2\theta`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> theta = np.linspace(0, np.pi, 2001)
    >>> round(directivity(theta, np.sin(theta) ** 2), 5)
    1.5
    """
    th = _vector("theta", theta, min_length=3)
    U = _vector("power_pattern", power_pattern, min_length=3)
    if th.shape != U.shape or np.any(np.diff(th) <= 0) or not np.isclose(th[0], 0) or not np.isclose(th[-1], np.pi):
        raise ValueError("theta must increase from 0 to pi, with one pattern value per angle")
    _nonnegative(power_pattern=U)
    if not np.any(U > 0):
        raise ValueError("power_pattern must not be identically zero")
    return float(2 * U.max() / trapezoid(U * np.sin(th), th))


def half_wave_dipole_pattern(theta):
    r"""Normalised far-field pattern of a thin half-wave dipole, :math:`F(\theta) = \cos(\tfrac\pi2\cos\theta)/\sin\theta`.

    Parameters
    ----------
    theta : float or array_like
        Angle from the dipole axis [rad].

    Returns
    -------
    float or ndarray
        Field-strength pattern [-], one broadside (:math:`\theta = \pi/2`)
        and zero along the axis; square it for the power pattern.

    Notes
    -----
    A sinusoidal current distribution is assumed. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/communications`.

    Examples
    --------
    >>> float(half_wave_dipole_pattern(np.pi / 2)), float(half_wave_dipole_pattern(0.0))
    (1.0, 0.0)
    """
    _finite(theta=theta)
    th = np.asarray(theta, dtype=float)
    s = np.sin(th)
    on_axis = np.abs(s) < 1e-12
    F = np.where(on_axis, 0.0, np.cos(np.pi / 2 * np.cos(th)) / np.where(on_axis, 1.0, s))
    return float(F) if F.ndim == 0 else F


def half_wave_dipole_directivity():
    r"""Directivity of a thin half-wave dipole, :math:`D = 4/\operatorname{Cin}(2\pi) \approx 1.641` (2.15 dBi).

    Returns
    -------
    float
        Directivity [-] (linear), from the closed-form integral of the
        pattern, :math:`\operatorname{Cin}(x) = \gamma + \ln x - \operatorname{Ci}(x)`.

    Notes
    -----
    See :doc:`/guides/communications`.

    Examples
    --------
    >>> round(half_wave_dipole_directivity(), 4)
    1.6409
    """
    return 4 / _cin_2pi()


def half_wave_dipole_resistance(impedance=FREE_SPACE_IMPEDANCE):
    r"""Radiation resistance of a thin half-wave dipole, :math:`R = \dfrac{\eta}{4\pi}\operatorname{Cin}(2\pi) \approx 73.1\ \Omega`.

    Parameters
    ----------
    impedance : float or array_like, optional
        Wave impedance of the medium :math:`\eta` [ohm], positive.

    Returns
    -------
    float or ndarray
        Radiation resistance referred to the feed-point current [ohm].

    Notes
    -----
    Carter's induced-EMF method with a sinusoidal current. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/communications`.

    References
    ----------
    P. S. Carter, "Circuit relations in radiating systems and applications
    to antenna problems," *Proc. IRE* 20 (1932) 1004-1041.

    Examples
    --------
    >>> round(half_wave_dipole_resistance(), 2)
    73.08
    """
    _positive(impedance=impedance)
    R = np.asarray(impedance, dtype=float) / (4 * np.pi) * _cin_2pi()
    return float(R) if R.ndim == 0 else R
