r"""Axial and torsional members: Hooke's law, Poisson's ratio and Coulomb's torsion.

A bar of length :math:`L`, area :math:`A` and modulus :math:`E` under axial
force :math:`P` stretches :math:`\delta = PL/(EA)` and stores
:math:`U = P^2L/(2EA)`. A circular shaft under torque :math:`T` carries
shear stress :math:`\tau = Tr/J` and twists :math:`\phi = TL/(GJ)`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "hooke_stress",
    "axial_deformation",
    "strain_energy_axial",
    "thermal_strain",
    "poisson_lateral_strain",
    "torsion_shear_stress",
    "torsion_angle",
]


def hooke_stress(E, strain):
    r"""Uniaxial stress from Hooke's law, :math:`\sigma = E\varepsilon`.

    Parameters
    ----------
    E : float or array_like
        Young's modulus [Pa], positive.
    strain : float or array_like
        Axial strain [-].

    Returns
    -------
    float or ndarray
        Stress [Pa], tension positive.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, strain=strain)
    _positive(E=E)
    return np.asarray(E) * np.asarray(strain)


def axial_deformation(force, length, area, E):
    r"""Elongation of a uniform elastic bar, :math:`\delta = PL/(AE)`.

    Parameters
    ----------
    force : float or array_like
        Axial force [N], tension positive.
    length : float or array_like
        Bar length [m], positive.
    area : float or array_like
        Cross-sectional area [m^2], positive.
    E : float or array_like
        Young's modulus [Pa], positive.

    Returns
    -------
    float or ndarray
        Elongation [m] (negative when shortened).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(force=force, length=length, area=area, E=E)
    _positive(length=length, area=area, E=E)
    return np.asarray(force) * length / (np.asarray(area) * E)


def strain_energy_axial(force, length, area, E):
    r"""Strain energy stored in a uniform elastic bar, :math:`U = P^2L/(2AE)`.

    Parameters
    ----------
    force : float or array_like
        Axial force [N].
    length : float or array_like
        Bar length [m], positive.
    area : float or array_like
        Cross-sectional area [m^2], positive.
    E : float or array_like
        Young's modulus [Pa], positive.

    Returns
    -------
    float or ndarray
        Strain energy [J], equal to :math:`P\delta/2`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(force=force, length=length, area=area, E=E)
    _positive(length=length, area=area, E=E)
    return np.asarray(force) ** 2 * length / (2 * np.asarray(area) * E)


def thermal_strain(alpha, delta_T):
    r"""Free thermal strain, :math:`\varepsilon_T = \alpha\,\Delta T`.

    Parameters
    ----------
    alpha : float or array_like
        Coefficient of linear thermal expansion [1/K].
    delta_T : float or array_like
        Temperature change [K].

    Returns
    -------
    float or ndarray
        Strain [-].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(alpha=alpha, delta_T=delta_T)
    return alpha * delta_T


def poisson_lateral_strain(axial_strain, nu):
    r"""Lateral strain of an isotropic bar, :math:`\varepsilon_t = -\nu\,\varepsilon_a`.

    Parameters
    ----------
    axial_strain : float or array_like
        Axial strain [-].
    nu : float or array_like
        Poisson's ratio [-], in the open interval (-1, 0.5) required for a
        positive-definite isotropic elastic material.

    Returns
    -------
    float or ndarray
        Lateral strain [-].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    S.-D. Poisson, "Mémoire sur l'équilibre et le mouvement des corps
    élastiques," *Mém. Acad. Roy. Sci. Inst. France* 8 (1829) 357-570.
    """
    _finite(axial_strain=axial_strain, nu=nu)
    _interval("nu", nu, -1, 0.5, closed=False)
    return -np.asarray(nu) * axial_strain


def torsion_shear_stress(T, r, J):
    r"""Shear stress in a circular shaft, :math:`\tau = Tr/J`.

    Parameters
    ----------
    T : float or array_like
        Torque [N m].
    r : float or array_like
        Radial position [m], nonnegative (the outer radius gives the peak).
    J : float or array_like
        Polar second moment of area [m^4], positive (:math:`\pi d^4/32` for
        a solid shaft).

    Returns
    -------
    float or ndarray
        Shear stress [Pa].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T=T, r=r, J=J)
    _nonnegative(r=r)
    _positive(J=J)
    return np.asarray(T) * r / J


def torsion_angle(T, L, G, J):
    r"""Angle of twist of a uniform circular shaft, :math:`\phi = TL/(GJ)`.

    Parameters
    ----------
    T : float or array_like
        Torque [N m].
    L : float or array_like
        Shaft length [m], positive.
    G : float or array_like
        Shear modulus [Pa], positive.
    J : float or array_like
        Polar second moment of area [m^4], positive.

    Returns
    -------
    float or ndarray
        Twist [rad].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    C.-A. Coulomb, "Recherches théoriques et expérimentales sur la force de
    torsion et sur l'élasticité des fils de métal," *Mém. Acad. Roy. Sci.*
    (1784, publ. 1787) 229-269.
    """
    _finite(T=T, L=L, G=G, J=J)
    _positive(L=L, G=G, J=J)
    return np.asarray(T) * L / (np.asarray(G) * J)
