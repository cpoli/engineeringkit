r"""Linear elasticity: stress, strain, the isotropic constants and thermal stress.

Young measured stiffness as a single modulus, :math:`\sigma = E\varepsilon`.
Lamé showed that an isotropic solid needs exactly two constants and that any
pair (:math:`E` and :math:`\nu`, :math:`\lambda` and :math:`\mu`, :math:`K`
and :math:`G`) gives all the others. Duhamel added the strain a temperature
change produces, :math:`\alpha\,\Delta T`, and the stress that appears when it
is prevented.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = [
    "normal_stress",
    "normal_strain",
    "hooke_stress",
    "shear_modulus",
    "bulk_modulus",
    "lame_parameters",
    "thermal_strain",
    "constrained_thermal_stress",
    "resilience_modulus",
    "toughness",
]


def _isotropic(E, nu):
    _positive(E=E)
    _interval("nu", nu, -1, 0.5, closed=False)


def normal_stress(force, area):
    r"""Engineering normal stress :math:`\sigma = F/A_0`.

    Parameters
    ----------
    force : float or array_like
        Axial force :math:`F` [N], positive in tension.
    area : float or array_like
        Original cross-sectional area :math:`A_0` [m^2], positive.

    Returns
    -------
    float or ndarray
        Normal stress [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> float(normal_stress(1000.0, 1e-4))
    10000000.0
    """
    _finite(force=force)
    _positive(area=area)
    return np.asarray(force) / area


def normal_strain(delta_length, length):
    r"""Engineering normal strain :math:`\varepsilon = \Delta L / L_0`.

    Parameters
    ----------
    delta_length : float or array_like
        Change of length :math:`\Delta L` [m], positive in extension.
    length : float or array_like
        Original gauge length :math:`L_0` [m], positive.

    Returns
    -------
    float or ndarray
        Strain [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _finite(delta_length=delta_length)
    _positive(length=length)
    return np.asarray(delta_length) / length


def hooke_stress(E, strain):
    r"""Uniaxial Hooke's law :math:`\sigma = E\varepsilon`.

    Parameters
    ----------
    E : float or array_like
        Young's modulus [Pa], positive.
    strain : float or array_like
        Normal strain [-].

    Returns
    -------
    float or ndarray
        Stress [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    T. Young, *A Course of Lectures on Natural Philosophy and the Mechanical
    Arts* (Johnson, London, 1807), vol. 1, lecture XIII.
    """
    _finite(strain=strain)
    _positive(E=E)
    return E * np.asarray(strain)


def shear_modulus(E, nu):
    r"""Isotropic shear modulus :math:`G = E / (2(1 + \nu))`.

    Parameters
    ----------
    E : float or array_like
        Young's modulus [Pa], positive.
    nu : float or array_like
        Poisson's ratio [-], in (-1, 0.5).

    Returns
    -------
    float or ndarray
        Shear modulus [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _isotropic(E, nu)
    return E / (2 * (1 + np.asarray(nu)))


def bulk_modulus(E, nu):
    r"""Isotropic bulk modulus :math:`K = E / (3(1 - 2\nu))`.

    Parameters
    ----------
    E : float or array_like
        Young's modulus [Pa], positive.
    nu : float or array_like
        Poisson's ratio [-], in (-1, 0.5); :math:`K \to \infty` as
        :math:`\nu \to 1/2`.

    Returns
    -------
    float or ndarray
        Bulk modulus [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _isotropic(E, nu)
    return E / (3 * (1 - 2 * np.asarray(nu)))


def lame_parameters(E, nu):
    r"""Lamé's constants :math:`\lambda` and :math:`\mu = G` from :math:`E` and :math:`\nu`.

    .. math:: \lambda = \frac{E\nu}{(1 + \nu)(1 - 2\nu)}, \qquad
              \mu = \frac{E}{2(1 + \nu)}

    so that :math:`\sigma_{ij} = \lambda\,\varepsilon_{kk}\delta_{ij} + 2\mu\,\varepsilon_{ij}`.

    Parameters
    ----------
    E : float or array_like
        Young's modulus [Pa], positive.
    nu : float or array_like
        Poisson's ratio [-], in (-1, 0.5).

    Returns
    -------
    tuple
        ``(lambda, mu)`` [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    G. Lamé, *Leçons sur la théorie mathématique de l'élasticité des corps
    solides* (Bachelier, Paris, 1852).
    """
    _isotropic(E, nu)
    nu = np.asarray(nu)
    return E * nu / ((1 + nu) * (1 - 2 * nu)), shear_modulus(E, nu)


def thermal_strain(alpha, dT):
    r"""Free thermal strain :math:`\varepsilon_T = \alpha\,\Delta T`.

    Parameters
    ----------
    alpha : float or array_like
        Linear expansion coefficient [1/K].
    dT : float or array_like
        Temperature change [K].

    Returns
    -------
    float or ndarray
        Strain [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _finite(alpha=alpha, dT=dT)
    return alpha * np.asarray(dT)


def constrained_thermal_stress(E, alpha, dT):
    r"""Stress in a bar held at fixed length while its temperature changes.

    .. math:: \sigma = -E\alpha\,\Delta T

    Parameters
    ----------
    E : float or array_like
        Young's modulus [Pa], positive.
    alpha : float or array_like
        Linear expansion coefficient [1/K].
    dT : float or array_like
        Temperature change [K].

    Returns
    -------
    float or ndarray
        Stress [Pa], compressive (negative) on heating.

    Notes
    -----
    Uniaxial restraint only; a plate restrained in two directions carries
    :math:`-E\alpha\Delta T/(1 - \nu)`. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    J.-M.-C. Duhamel, "Mémoire sur le calcul des actions moléculaires
    développées par les changements de température dans les corps solides,"
    *Mém. Divers Savants Acad. R. Sci. Inst. France* 5 (1838) 440-498.
    """
    _finite(alpha=alpha, dT=dT)
    _positive(E=E)
    return -E * alpha * np.asarray(dT)


def resilience_modulus(sigma_y, E):
    r"""Elastic strain energy per volume at yield, :math:`U_r = \sigma_y^2 / 2E`.

    Parameters
    ----------
    sigma_y : float or array_like
        Yield strength [Pa], positive.
    E : float or array_like
        Young's modulus [Pa], positive.

    Returns
    -------
    float or ndarray
        Modulus of resilience [J/m^3].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _positive(sigma_y=sigma_y, E=E)
    return np.asarray(sigma_y) ** 2 / (2 * E)


def toughness(strain, stress):
    r"""Area under a stress-strain curve, :math:`\int \sigma\,d\varepsilon`.

    Parameters
    ----------
    strain : array_like
        Strictly increasing strains [-], at least two.
    stress : array_like
        Matching stresses [Pa].

    Returns
    -------
    float
        Energy absorbed per unit volume [J/m^3], by the trapezoid rule.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    strain, stress = np.asarray(strain, float), np.asarray(stress, float)
    _finite(strain=strain, stress=stress)
    if strain.ndim != 1 or strain.shape != stress.shape or len(strain) < 2 or np.any(np.diff(strain) <= 0):
        raise ValueError("stress and increasing strain must be matching vectors")
    return float(trapezoid(stress, strain))
