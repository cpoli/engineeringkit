r"""Shafts under combined bending and torsion.

A solid round shaft carrying bending moment :math:`M` and torque :math:`T`
has, at its surface, the bending stress :math:`\sigma = 32M/\pi d^3` and the
torsional shear :math:`\tau = 16T/\pi d^3`. The distortion-energy (von
Mises) and maximum-shear (Tresca, Guest) criteria combine them into one
equivalent stress, and under fatigue the Goodman line weighs the
alternating part against the endurance limit and the mean part against the
ultimate strength.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = ["shaft_von_mises_stress", "shaft_max_shear_stress", "shaft_diameter", "shaft_diameter_goodman"]

_CRITERIA = ("von_mises", "tresca")


def shaft_von_mises_stress(d, M, T):
    r"""Von Mises stress at the surface of a solid shaft, :math:`\sigma' = \frac{16}{\pi d^3}\sqrt{4M^2 + 3T^2}`.

    Parameters
    ----------
    d : float or array_like
        Shaft diameter [m], positive.
    M : float or array_like
        Bending moment [N m].
    T : float or array_like
        Torque [N m].

    Returns
    -------
    float or ndarray
        :math:`\sqrt{\sigma^2 + 3\tau^2}` [Pa].

    Notes
    -----
    Transverse shear, which is zero at the surface fibre of greatest bending
    stress, and axial load are neglected. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(shaft_von_mises_stress(0.04, 500.0, 0.0)) / 1e6, 4)
    79.5775
    """
    _finite(M=M, T=T)
    _positive(d=d)
    return 16 / (np.pi * np.asarray(d, dtype=float) ** 3) * np.sqrt(4 * np.square(M) + 3 * np.square(T))


def shaft_max_shear_stress(d, M, T):
    r"""Maximum shear stress at the surface of a solid shaft, :math:`\tau_\text{max} = \frac{16}{\pi d^3}\sqrt{M^2 + T^2}`.

    Parameters
    ----------
    d : float or array_like
        Shaft diameter [m], positive.
    M : float or array_like
        Bending moment [N m].
    T : float or array_like
        Torque [N m].

    Returns
    -------
    float or ndarray
        :math:`\sqrt{(\sigma/2)^2 + \tau^2}` [Pa], the radius of Mohr's circle.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    J. J. Guest, "On the strength of ductile materials under combined
    stress," *Phil. Mag.* (5) 50 (1900) 69-132.

    Examples
    --------
    >>> round(float(shaft_max_shear_stress(0.04, 300.0, 400.0)) / 1e6, 4)
    39.7887
    """
    _finite(M=M, T=T)
    _positive(d=d)
    return 16 / (np.pi * np.asarray(d, dtype=float) ** 3) * np.hypot(M, T)


def shaft_diameter(M, T, yield_strength, safety_factor=1.0, criterion="von_mises"):
    r"""Smallest solid-shaft diameter that keeps static bending plus torsion below yield.

    .. math::

        d = \left(\frac{16 n}{\pi S_y}\sqrt{4M^2 + 3T^2}\right)^{1/3}
        \ \text{(von Mises)}, \qquad
        d = \left(\frac{32 n}{\pi S_y}\sqrt{M^2 + T^2}\right)^{1/3}
        \ \text{(Tresca)}

    Parameters
    ----------
    M : float or array_like
        Bending moment [N m].
    T : float or array_like
        Torque [N m].
    yield_strength : float or array_like
        Tensile yield strength :math:`S_y` [Pa], positive.
    safety_factor : float or array_like, optional
        :math:`n` [-], positive.
    criterion : {"von_mises", "tresca"}, optional
        Yield criterion; Tresca is the more conservative.

    Returns
    -------
    float or ndarray
        Diameter [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(shaft_diameter(500.0, 0.0, 79.5775e6)), 6)
    0.04
    """
    _finite(M=M, T=T)
    _positive(yield_strength=yield_strength, safety_factor=safety_factor)
    if criterion not in _CRITERIA:
        raise ValueError(f"criterion must be one of {_CRITERIA}")
    if criterion == "von_mises":
        demand = 16 * np.sqrt(4 * np.square(M) + 3 * np.square(T))
    else:
        demand = 32 * np.hypot(M, T)
    return np.cbrt(safety_factor * demand / (np.pi * np.asarray(yield_strength, dtype=float)))


def shaft_diameter_goodman(endurance_limit, ultimate_strength, Ma=0.0, Tm=0.0, Mm=0.0, Ta=0.0, safety_factor=1.0, Kf=1.0, Kfs=1.0):
    r"""Shaft diameter for infinite fatigue life by the distortion-energy Goodman criterion.

    .. math::

        d = \left\{\frac{16 n}{\pi}\left[
        \frac{\sqrt{4(K_f M_a)^2 + 3(K_{fs} T_a)^2}}{S_e} +
        \frac{\sqrt{4(K_f M_m)^2 + 3(K_{fs} T_m)^2}}{S_{ut}}\right]\right\}^{1/3}

    Parameters
    ----------
    endurance_limit : float
        Fully corrected endurance limit :math:`S_e` [Pa], positive.
    ultimate_strength : float
        Ultimate tensile strength :math:`S_{ut}` [Pa], positive.
    Ma, Mm : float, optional
        Alternating and mean bending moment [N m]. A rotating shaft under a
        steady bending load has only :math:`M_a`.
    Ta, Tm : float, optional
        Alternating and mean torque [N m].
    safety_factor : float, optional
        :math:`n` [-], positive.
    Kf, Kfs : float, optional
        Fatigue stress-concentration factors in bending and torsion [-], at
        least one.

    Returns
    -------
    float
        Diameter [m] at which the von Mises alternating and mean stresses lie
        on the Goodman line :math:`\sigma_a'/S_e + \sigma_m'/S_{ut} = 1/n`.

    Notes
    -----
    Yielding on the first cycle is not checked; compare
    :math:`\sigma_a' + \sigma_m'` with :math:`S_y` separately. Inputs are raw
    SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(shaft_diameter_goodman(200e6, 600e6, Ma=100.0, Tm=0.0), 6)
    0.017205
    """
    _positive(endurance_limit=endurance_limit, ultimate_strength=ultimate_strength, safety_factor=safety_factor)
    _finite(Ma=Ma, Tm=Tm, Mm=Mm, Ta=Ta)
    _interval("Kf", Kf, 1, np.inf, closed="left")
    _interval("Kfs", Kfs, 1, np.inf, closed="left")
    alternating = np.sqrt(4 * (Kf * Ma) ** 2 + 3 * (Kfs * Ta) ** 2) / endurance_limit
    mean = np.sqrt(4 * (Kf * Mm) ** 2 + 3 * (Kfs * Tm) ** 2) / ultimate_strength
    if alternating + mean == 0:
        raise ValueError("at least one load must be nonzero")
    return float(np.cbrt(16 * safety_factor / np.pi * (alternating + mean)))
