r"""Two-dimensional supersonic flow: oblique shocks and Prandtl-Meyer expansions.

A supersonic stream turned into itself by a wedge of half-angle
:math:`\theta` passes through a straight oblique shock at the wave angle
:math:`\beta` given by the :math:`\theta`-:math:`\beta`-:math:`M` relation;
the jump is a normal shock in the velocity component normal to the wave. A
stream turned away from itself expands isentropically through a centred
fan, and the Prandtl-Meyer function :math:`\nu(M)` gives the turning angle.
Both were worked out by Theodor Meyer in his 1908 Göttingen thesis.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.aerodynamics.core.base import ObliqueShockResult
from engineeringkit.aerodynamics.systems.compressible import (
    mach_angle,
    normal_shock_mach2,
    normal_shock_pressure_ratio,
    stagnation_pressure,
)
from engineeringkit.aerodynamics.utils.gas import check_gamma

__all__ = [
    "oblique_shock_deflection",
    "oblique_shock_max_deflection",
    "oblique_shock_angle",
    "oblique_shock",
    "prandtl_meyer_angle",
    "inverse_prandtl_meyer",
    "prandtl_meyer_expansion",
]


def oblique_shock_deflection(M1, beta, gamma=1.4):
    r"""Flow deflection behind an oblique shock, the :math:`\theta`-:math:`\beta`-:math:`M` relation.

    .. math::

        \tan\theta = 2\cot\beta\,\frac{M_1^2\sin^2\beta - 1}{M_1^2(\gamma + \cos 2\beta) + 2}

    Parameters
    ----------
    M1 : float or array_like
        Upstream Mach number [-], at least one.
    beta : float or array_like
        Wave angle from the upstream flow [rad], in
        :math:`[\mu, \pi/2]` where :math:`\mu = \arcsin(1/M_1)`.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Deflection angle [rad]; zero for a Mach wave (:math:`\beta = \mu`)
        and for a normal shock (:math:`\beta = \pi/2`).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Modern Compressible Flow*, 3rd ed. (McGraw-Hill,
    2003), eq. 4.17.

    Examples
    --------
    >>> round(float(np.degrees(oblique_shock_deflection(2.0, np.radians(53.4229)))), 3)
    20.0
    """
    _interval("M1", M1, 1, np.inf, closed="left")
    check_gamma(gamma)
    M1, beta = np.asarray(M1, float), np.asarray(beta, float)
    _interval("beta", beta, 0, np.pi / 2, closed="right")
    if np.any(beta < mach_angle(M1) - 1e-12):
        raise ValueError("beta must be at least the Mach angle arcsin(1/M1)")
    mn1sq = np.maximum(M1**2 * np.sin(beta) ** 2 - 1, 0.0)
    return np.arctan2(2 * mn1sq * np.cos(beta), np.sin(beta) * (M1**2 * (gamma + np.cos(2 * beta)) + 2))


def oblique_shock_max_deflection(M1, gamma=1.4):
    r"""Largest deflection an attached oblique shock can turn, and its wave angle.

    The wave angle at maximum deflection satisfies

    .. math::

        \sin^2\beta_{\max} = \frac{1}{\gamma M_1^2}\left[\frac{\gamma+1}{4}M_1^2 - 1
        + \sqrt{(\gamma+1)\left(\frac{\gamma+1}{16}M_1^4 + \frac{\gamma-1}{2}M_1^2 + 1\right)}\right].

    Parameters
    ----------
    M1 : float or array_like
        Upstream Mach number [-], at least one.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    theta_max : float or ndarray
        Maximum deflection [rad]; a steeper wedge carries a detached bow shock.
    beta_max : float or ndarray
        Wave angle at that deflection [rad], dividing the weak and strong branches.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> theta_max, beta_max = oblique_shock_max_deflection(2.0)
    >>> round(float(np.degrees(theta_max)), 2), round(float(np.degrees(beta_max)), 2)
    (22.97, 64.67)
    """
    _interval("M1", M1, 1, np.inf, closed="left")
    check_gamma(gamma)
    m2 = np.asarray(M1, float) ** 2
    root = np.sqrt((gamma + 1) * ((gamma + 1) / 16 * m2**2 + (gamma - 1) / 2 * m2 + 1))
    beta_max = np.arcsin(np.sqrt(((gamma + 1) / 4 * m2 - 1 + root) / (gamma * m2)))
    return oblique_shock_deflection(M1, beta_max, gamma), beta_max


def oblique_shock_angle(M1, theta, gamma=1.4, strong=False):
    r"""Wave angle :math:`\beta` of the attached oblique shock that turns a stream by :math:`\theta`.

    Solves the :math:`\theta`-:math:`\beta`-:math:`M` relation, a cubic in
    :math:`\tan\beta`, in closed form:

    .. math::

        \tan\beta = \frac{M_1^2 - 1 + 2\lambda\cos\left[(4\pi\delta + \cos^{-1}\chi)/3\right]}
        {3\left(1 + \tfrac{\gamma-1}{2}M_1^2\right)\tan\theta},

    with :math:`\delta = 1` for the weak and :math:`\delta = 0` for the strong
    solution.

    Parameters
    ----------
    M1 : float or array_like
        Upstream Mach number [-], greater than one.
    theta : float or array_like
        Deflection angle [rad], positive and at most
        :func:`oblique_shock_max_deflection`.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.
    strong : bool, optional
        Return the strong-shock branch (subsonic downstream flow) instead of
        the weak one that occurs on wedges and ramps.

    Returns
    -------
    float or ndarray
        Wave angle [rad].

    Raises
    ------
    ValueError
        If ``theta`` exceeds the maximum deflection (the shock detaches).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    M. J. Thompson, "A note on the calculation of oblique shock-wave
    characteristics," *J. Aeronaut. Sci.* 17 (1950) 744; T. Meyer, "Über
    zweidimensionale Bewegungsvorgänge in einem Gas, das mit
    Überschallgeschwindigkeit strömt," *Forschungsheft VDI* 62 (1908).

    Examples
    --------
    >>> round(float(np.degrees(oblique_shock_angle(2.0, np.radians(20.0)))), 2)
    53.42
    """
    _interval("M1", M1, 1, np.inf, closed=False)
    _positive(theta=theta)
    check_gamma(gamma)
    M1, theta = np.asarray(M1, float), np.asarray(theta, float)
    if np.any(theta > oblique_shock_max_deflection(M1, gamma)[0] * (1 + 1e-12)):
        raise ValueError("theta exceeds the maximum deflection: the shock is detached")
    m2, tan2 = M1**2, np.tan(theta) ** 2
    k1 = 1 + (gamma - 1) / 2 * m2
    lam = np.sqrt(np.maximum((m2 - 1) ** 2 - 3 * k1 * (1 + (gamma + 1) / 2 * m2) * tan2, 0.0))
    chi = ((m2 - 1) ** 3 - 9 * k1 * (k1 + (gamma + 1) / 4 * m2**2) * tan2) / lam**3
    delta = 0 if strong else 1
    tan_beta = (m2 - 1 + 2 * lam * np.cos((4 * np.pi * delta + np.arccos(np.clip(chi, -1, 1))) / 3)) / (3 * k1 * np.tan(theta))
    return np.arctan(tan_beta)


def oblique_shock(M1, theta, gamma=1.4, strong=False) -> ObliqueShockResult:
    r"""Downstream state of an attached oblique shock.

    The velocity component normal to the wave, :math:`M_{n1} = M_1\sin\beta`,
    passes through a normal shock; the tangential component is unchanged, so
    :math:`M_2 = M_{n2}/\sin(\beta - \theta)`.

    Parameters
    ----------
    M1 : float
        Upstream Mach number [-], greater than one.
    theta : float
        Deflection angle [rad], positive and at most the maximum deflection.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.
    strong : bool, optional
        Use the strong-shock branch instead of the weak one.

    Returns
    -------
    ObliqueShockResult
        Wave angle, downstream Mach number and the static and stagnation
        ratios across the shock.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> shock = oblique_shock(2.0, np.radians(20.0))
    >>> round(shock.mach2, 3), round(shock.pressure_ratio, 3)
    (1.21, 2.843)
    """
    beta = float(oblique_shock_angle(M1, theta, gamma, strong))
    mn1 = M1 * np.sin(beta)
    mn2 = float(normal_shock_mach2(mn1, gamma))
    p_ratio = float(normal_shock_pressure_ratio(mn1, gamma))
    rho_ratio = (gamma + 1) * mn1**2 / ((gamma - 1) * mn1**2 + 2)
    mach2 = mn2 / np.sin(beta - theta)
    p0_ratio = float(stagnation_pressure(p_ratio, mach2, gamma) / stagnation_pressure(1.0, M1, gamma))
    return ObliqueShockResult(
        mach1=float(M1),
        theta=float(theta),
        beta=beta,
        mach2=float(mach2),
        pressure_ratio=p_ratio,
        density_ratio=float(rho_ratio),
        temperature_ratio=float(p_ratio / rho_ratio),
        stagnation_pressure_ratio=p0_ratio,
        gamma=float(gamma),
        strong=bool(strong),
    )


def prandtl_meyer_angle(M, gamma=1.4):
    r"""Prandtl-Meyer function :math:`\nu(M)`, the turning angle that expands a sonic stream to :math:`M`.

    .. math::

        \nu(M) = \sqrt{\frac{\gamma+1}{\gamma-1}}\arctan\sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}
        - \arctan\sqrt{M^2-1}

    Parameters
    ----------
    M : float or array_like
        Mach number [-], at least one.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Prandtl-Meyer angle [rad]; zero at :math:`M = 1`, tending to
        :math:`\tfrac{\pi}{2}(\sqrt{(\gamma+1)/(\gamma-1)} - 1)` as
        :math:`M \to \infty`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    L. Prandtl, "Neue Untersuchungen über die strömende Bewegung der Gase
    und Dämpfe," *Physikalische Zeitschrift* 8 (1907) 23-30; T. Meyer,
    *Forschungsheft VDI* 62 (1908).

    Examples
    --------
    >>> round(float(np.degrees(prandtl_meyer_angle(2.0))), 3)
    26.38
    """
    _interval("M", M, 1, np.inf, closed="left")
    check_gamma(gamma)
    root = np.sqrt(np.asarray(M, float) ** 2 - 1)
    k = np.sqrt((gamma + 1) / (gamma - 1))
    return k * np.arctan(root / k) - np.arctan(root)


def inverse_prandtl_meyer(nu, gamma=1.4):
    r"""Mach number whose Prandtl-Meyer angle is :math:`\nu`, found by a bracketing root search.

    Parameters
    ----------
    nu : float or array_like
        Prandtl-Meyer angle [rad], nonnegative and below the limit
        :math:`\tfrac{\pi}{2}(\sqrt{(\gamma+1)/(\gamma-1)} - 1)` (130.45° for
        :math:`\gamma = 1.4`).
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Mach number [-], at least one.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(inverse_prandtl_meyer(np.radians(26.38))), 4)
    2.0
    """
    check_gamma(gamma)
    nu_max = np.pi / 2 * (np.sqrt((gamma + 1) / (gamma - 1)) - 1)
    _interval("nu", nu, 0, nu_max, closed="left")

    def solve(target):
        if target == 0:
            return 1.0
        upper = 2.0
        while prandtl_meyer_angle(upper, gamma) < target:
            upper *= 2
        return brentq(lambda m: prandtl_meyer_angle(m, gamma) - target, 1.0, upper, xtol=1e-14, rtol=1e-14)

    result = np.vectorize(solve, otypes=[float])(np.asarray(nu, float))
    return result[()] if result.ndim == 0 else result


def prandtl_meyer_expansion(M1, theta, gamma=1.4):
    r"""Mach number after a supersonic stream turns away from itself by :math:`\theta`.

    The expansion is isentropic, so :math:`\nu(M_2) = \nu(M_1) + \theta`
    and the stagnation state is unchanged.

    Parameters
    ----------
    M1 : float or array_like
        Upstream Mach number [-], at least one.
    theta : float or array_like
        Turning angle [rad], nonnegative.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Downstream Mach number [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(prandtl_meyer_expansion(1.5, np.radians(20.0))), 3)
    2.207
    """
    _nonnegative(theta=theta)
    return inverse_prandtl_meyer(prandtl_meyer_angle(M1, gamma) + np.asarray(theta, float), gamma)
