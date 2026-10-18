r"""Fracture mechanics: Griffith's energy balance, Irwin's K and Paris's law.

Griffith showed that a crack grows when the elastic energy it releases pays
for the new surface, which makes strength fall as :math:`1/\sqrt{a}`. Irwin
recast the balance as the stress intensity factor
:math:`K = Y\sigma\sqrt{\pi a}` reaching the toughness :math:`K_{Ic}`. Paris
found that a fatigue crack grows by a power of the range of :math:`K` each
cycle.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.integrate import cumulative_trapezoid

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.materials.core.base import CrackGrowthResult

__all__ = ["griffith_fracture_stress", "stress_intensity", "critical_crack_length", "paris_crack_growth", "crack_growth_life"]


def griffith_fracture_stress(E, gamma, a):
    r"""Griffith's fracture stress of a brittle plate with a central crack.

    .. math:: \sigma_f = \sqrt{\frac{2E\gamma_s}{\pi a}}

    Parameters
    ----------
    E : float or array_like
        Young's modulus [Pa], positive.
    gamma : float or array_like
        Surface energy :math:`\gamma_s` [J/m^2], positive.
    a : float or array_like
        Half-length of the crack [m], positive.

    Returns
    -------
    float or ndarray
        Remote tensile stress at fracture [Pa], plane stress.

    Notes
    -----
    Ideally brittle solid; plastic work at the tip, which dominates in
    metals, is left out. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/materials`.

    References
    ----------
    A. A. Griffith, "The phenomena of rupture and flow in solids," *Phil.
    Trans. R. Soc. Lond. A* 221 (1921) 163-198.
    """
    _positive(E=E, gamma=gamma, a=a)
    return np.sqrt(2 * E * gamma / (np.pi * np.asarray(a)))


def stress_intensity(K_geometry, sigma, a):
    r"""Mode-I stress intensity factor :math:`K_I = Y\sigma\sqrt{\pi a}`.

    Parameters
    ----------
    K_geometry : float or array_like
        Geometry factor :math:`Y` [-], positive (1 for a central crack in an
        infinite plate, about 1.12 for an edge crack).
    sigma : float or array_like
        Remote tensile stress [Pa], nonnegative.
    a : float or array_like
        Crack length (edge) or half-length (central) [m], positive.

    Returns
    -------
    float or ndarray
        Stress intensity [Pa sqrt(m)].

    Notes
    -----
    Linear-elastic fracture mechanics, valid while the plastic zone is small
    compared with the crack. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/materials`.

    References
    ----------
    G. R. Irwin, "Analysis of stresses and strains near the end of a crack
    traversing a plate," *J. Appl. Mech.* 24 (1957) 361-364.
    """
    _positive(K_geometry=K_geometry, a=a)
    _nonnegative(sigma=sigma)
    return K_geometry * np.asarray(sigma) * np.sqrt(np.pi * np.asarray(a))


def critical_crack_length(Kic, sigma, Y=1.0):
    r"""Crack length at which :math:`K_I` reaches the toughness, :math:`a_c = (K_{Ic}/Y\sigma)^2/\pi`.

    Parameters
    ----------
    Kic : float or array_like
        Plane-strain fracture toughness [Pa sqrt(m)], positive.
    sigma : float or array_like
        Remote tensile stress [Pa], positive.
    Y : float or array_like
        Geometry factor [-], positive.

    Returns
    -------
    float or ndarray
        Critical crack length [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _positive(Kic=Kic, sigma=sigma, Y=Y)
    return (Kic / (Y * np.asarray(sigma))) ** 2 / np.pi


def paris_crack_growth(C, m, delta_K):
    r"""Fatigue crack growth rate from Paris's law, :math:`da/dN = C\,\Delta K^m`.

    Parameters
    ----------
    C : float or array_like
        Paris coefficient [m/cycle per (Pa sqrt(m))^m], nonnegative. Values
        quoted for MPa sqrt(m) must be multiplied by :math:`10^{-6m}`.
    m : float or array_like
        Paris exponent [-], positive (about 3 for steels).
    delta_K : float or array_like
        Stress intensity range [Pa sqrt(m)], nonnegative.

    Returns
    -------
    float or ndarray
        Growth per cycle [m/cycle].

    Notes
    -----
    The stable (Region II) regime only: neither the threshold
    :math:`\Delta K_{th}` nor the acceleration near :math:`K_{Ic}` is
    modelled. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    P. Paris and F. Erdogan, "A critical analysis of crack propagation laws,"
    *J. Basic Eng.* 85 (1963) 528-533.
    """
    _nonnegative(C=C, delta_K=delta_K)
    _positive(m=m)
    return C * np.asarray(delta_K) ** m


def crack_growth_life(C, m, delta_sigma, a0, af, Y: float | Callable = 1.0, points=2001):
    r"""Integrate Paris's law from an initial to a final crack length.

    .. math:: N = \int_{a_0}^{a_f} \frac{da}{C\,(Y(a)\,\Delta\sigma\sqrt{\pi a})^m}

    Parameters
    ----------
    C : float
        Paris coefficient [m/cycle per (Pa sqrt(m))^m], positive.
    m : float
        Paris exponent [-], positive.
    delta_sigma : float
        Stress range [Pa], positive.
    a0 : float
        Initial crack length [m], positive.
    af : float
        Final crack length [m], greater than ``a0``; often
        :func:`critical_crack_length` at the peak stress.
    Y : float or callable
        Geometry factor [-], a positive constant or a function ``Y(a)`` of
        crack length [m] (for example a finite-width correction).
    points : int
        Number of geometrically spaced crack lengths, at least 2.

    Returns
    -------
    CrackGrowthResult
        Crack length against cycles.

    Notes
    -----
    The trapezoid rule on a geometric grid; with a constant ``Y`` and
    :math:`m \ne 2` the closed form
    :math:`N = 2(a_0^{1-m/2} - a_f^{1-m/2}) / ((m - 2) C (Y\Delta\sigma\sqrt\pi)^m)`
    checks it. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _positive(C=C, m=m, delta_sigma=delta_sigma, a0=a0, af=af)
    _integer("points", points, 2)
    if af <= a0:
        raise ValueError("af must exceed a0")
    a = np.geomspace(a0, af, points)
    geometry = np.asarray(Y(a), float) * np.ones_like(a) if callable(Y) else np.full_like(a, Y)
    _positive(Y=geometry)
    rate = paris_crack_growth(C, m, stress_intensity(geometry, delta_sigma, a))
    cycles = cumulative_trapezoid(1 / rate, a, initial=0)
    return CrackGrowthResult(crack_length=a, cycles=cycles)
