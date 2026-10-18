r"""Rules of mixtures: the Voigt and Reuss bounds on a composite's stiffness.

Voigt assumed every phase strains equally (fibres loaded along their length),
which averages the moduli, :math:`E_V = V_f E_f + (1 - V_f) E_m`. Reuss assumed
every phase carries the same stress (loaded across), which averages the
compliances. Hill (1952) proved that the true modulus of any arrangement of
the phases lies between the two.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = ["rule_of_mixtures_E_longitudinal", "rule_of_mixtures_E_transverse", "density_rule_of_mixtures"]


def rule_of_mixtures_E_longitudinal(Ef, Em, Vf):
    r"""Voigt (equal-strain) modulus :math:`E_1 = V_f E_f + (1 - V_f) E_m`.

    Parameters
    ----------
    Ef : float or array_like
        Fibre modulus [Pa], positive.
    Em : float or array_like
        Matrix modulus [Pa], positive.
    Vf : float or array_like
        Fibre volume fraction [-], in [0, 1].

    Returns
    -------
    float or ndarray
        Modulus along continuous aligned fibres [Pa], the upper bound.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    W. Voigt, "Ueber die Beziehung zwischen den beiden Elasticitätsconstanten
    isotroper Körper," *Ann. Phys.* 274 (1889) 573-587.
    """
    _positive(Ef=Ef, Em=Em)
    _interval("Vf", Vf, 0, 1)
    Vf = np.asarray(Vf)
    return Ef * Vf + Em * (1 - Vf)


def rule_of_mixtures_E_transverse(Ef, Em, Vf):
    r"""Reuss (equal-stress) modulus :math:`1/E_2 = V_f/E_f + (1 - V_f)/E_m`.

    Parameters
    ----------
    Ef : float or array_like
        Fibre modulus [Pa], positive.
    Em : float or array_like
        Matrix modulus [Pa], positive.
    Vf : float or array_like
        Fibre volume fraction [-], in [0, 1].

    Returns
    -------
    float or ndarray
        Modulus across the fibres [Pa], the lower bound.

    Notes
    -----
    Real transverse moduli lie above the Reuss value because the matrix
    between neighbouring fibres is not uniformly stressed. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    A. Reuss, "Berechnung der Fließgrenze von Mischkristallen auf Grund der
    Plastizitätsbedingung für Einkristalle," *Z. Angew. Math. Mech.* 9 (1929)
    49-58.
    """
    _positive(Ef=Ef, Em=Em)
    _interval("Vf", Vf, 0, 1)
    Vf = np.asarray(Vf)
    return 1 / (Vf / Ef + (1 - Vf) / Em)


def density_rule_of_mixtures(rhof, rhom, Vf):
    r"""Composite density :math:`\rho = V_f\rho_f + (1 - V_f)\rho_m`.

    Parameters
    ----------
    rhof : float or array_like
        Fibre density [kg/m^3], positive.
    rhom : float or array_like
        Matrix density [kg/m^3], positive.
    Vf : float or array_like
        Fibre volume fraction [-], in [0, 1].

    Returns
    -------
    float or ndarray
        Density [kg/m^3], exact for a void-free composite.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _positive(rhof=rhof, rhom=rhom)
    _interval("Vf", Vf, 0, 1)
    Vf = np.asarray(Vf)
    return rhof * Vf + rhom * (1 - Vf)
