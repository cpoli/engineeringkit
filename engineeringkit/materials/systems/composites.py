r"""Composites: rules of mixtures, Halpin-Tsai and classical laminate theory.

Voigt assumed every phase strains equally (fibres loaded along their length),
which averages the moduli, :math:`E_V = V_f E_f + (1 - V_f) E_m`. Reuss assumed
every phase carries the same stress (loaded across), which averages the
compliances. Hill (1952) proved that the true modulus of any arrangement of
the phases lies between the two, and Halpin and Tsai (1969) interpolated
between them with one geometric parameter. Classical laminate theory stacks
unidirectional plies at different angles and relates in-plane forces and
moments to midplane strains and curvatures through the ABD matrix.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.materials.core.base import Laminate

__all__ = [
    "rule_of_mixtures_E_longitudinal",
    "rule_of_mixtures_E_transverse",
    "density_rule_of_mixtures",
    "halpin_tsai",
    "lamina_stiffness",
    "transformed_stiffness",
    "classical_laminate",
]


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


def halpin_tsai(Pf, Pm, Vf, xi):
    r"""Halpin-Tsai estimate of a composite modulus.

    .. math:: \frac{P}{P_m} = \frac{1 + \xi\eta V_f}{1 - \eta V_f}, \qquad
              \eta = \frac{P_f/P_m - 1}{P_f/P_m + \xi}

    Parameters
    ----------
    Pf : float or array_like
        Fibre modulus [Pa], positive (Young's or shear).
    Pm : float or array_like
        Matrix modulus [Pa], positive (the same kind).
    Vf : float or array_like
        Fibre volume fraction [-], in [0, 1].
    xi : float or array_like
        Reinforcement parameter [-], nonnegative: 2 for the transverse
        modulus :math:`E_2` of circular fibres, 1 for :math:`G_{12}`,
        :math:`2\ell/d` for :math:`E_1` of short fibres.

    Returns
    -------
    float or ndarray
        Composite modulus [Pa]. :math:`\xi = 0` gives the Reuss bound and
        :math:`\xi \to \infty` the Voigt bound.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    J. C. Halpin and S. W. Tsai, "Effects of environmental factors on
    composite materials," AFML-TR-67-423, Air Force Materials Laboratory
    (1969); J. C. Halpin and J. L. Kardos, "The Halpin-Tsai equations: a
    review," *Polym. Eng. Sci.* 16 (1976) 344-352.

    Examples
    --------
    >>> round(float(halpin_tsai(72e9, 3.5e9, 0.6, 2.0)) / 1e9, 2)
    14.89
    """
    _positive(Pf=Pf, Pm=Pm)
    _interval("Vf", Vf, 0, 1)
    _nonnegative(xi=xi)
    ratio, Vf = np.asarray(Pf) / Pm, np.asarray(Vf)
    eta = (ratio - 1) / (ratio + xi)
    return Pm * (1 + xi * eta * Vf) / (1 - eta * Vf)


def lamina_stiffness(E1, E2, nu12, G12):
    r"""Reduced stiffness matrix :math:`\mathbf Q` of an orthotropic ply in plane stress.

    .. math:: Q_{11} = \frac{E_1}{1 - \nu_{12}\nu_{21}}, \quad
              Q_{22} = \frac{E_2}{1 - \nu_{12}\nu_{21}}, \quad
              Q_{12} = \nu_{12} Q_{22}, \quad Q_{66} = G_{12}

    with :math:`\nu_{21} = \nu_{12} E_2/E_1`.

    Parameters
    ----------
    E1, E2 : float
        Moduli along and across the fibres [Pa], positive.
    nu12 : float
        Major Poisson's ratio [-], with :math:`\nu_{12}^2 < E_1/E_2`.
    G12 : float
        In-plane shear modulus [Pa], positive.

    Returns
    -------
    ndarray
        :math:`3 \times 3` matrix [Pa] relating
        :math:`(\sigma_1, \sigma_2, \tau_{12})` to
        :math:`(\varepsilon_1, \varepsilon_2, \gamma_{12})`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _positive(E1=E1, E2=E2, G12=G12)
    _finite(nu12=nu12)
    nu21 = nu12 * E2 / E1
    if nu12 * nu21 >= 1:
        raise ValueError("nu12**2 must be smaller than E1/E2 for a stable ply")
    d = 1 - nu12 * nu21
    return np.array([[E1 / d, nu12 * E2 / d, 0.0], [nu12 * E2 / d, E2 / d, 0.0], [0.0, 0.0, G12]], dtype=float)


def transformed_stiffness(Q, theta):
    r"""Ply stiffness :math:`\bar{\mathbf Q}` in laminate axes, for fibres at angle :math:`\theta` to x.

    .. math:: \bar{\mathbf Q} = \mathbf T^{-1}\mathbf Q\,\mathbf R\,\mathbf T\,\mathbf R^{-1}

    with :math:`\mathbf T` the stress transformation and :math:`\mathbf R`
    the Reuter matrix that converts engineering to tensor shear strain.

    Parameters
    ----------
    Q : array_like
        :math:`3 \times 3` reduced stiffness [Pa], from :func:`lamina_stiffness`.
    theta : float
        Fibre angle from the laminate x axis [rad], counterclockwise.

    Returns
    -------
    ndarray
        :math:`3 \times 3` matrix [Pa] relating
        :math:`(\sigma_x, \sigma_y, \tau_{xy})` to
        :math:`(\varepsilon_x, \varepsilon_y, \gamma_{xy})`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    Q = np.asarray(Q, float)
    if Q.shape != (3, 3):
        raise ValueError("Q must be a 3x3 matrix")
    _finite(Q=Q, theta=theta)
    T = _stress_rotation(theta)
    R = np.diag([1.0, 1.0, 2.0])
    return np.linalg.solve(T, Q @ R @ T @ np.linalg.inv(R))


def _stress_rotation(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c * c, s * s, 2 * s * c], [s * s, c * c, -2 * s * c], [-s * c, s * c, c * c - s * s]])


def classical_laminate(E1, E2, nu12, G12, angles, ply_thickness):
    r"""Assemble the ABD matrix of a laminate by classical laminate theory.

    .. math:: (A, B, D)_{ij} = \sum_k \bar Q_{ij}^{(k)}\,
              \left(z_k - z_{k-1},\ \tfrac12(z_k^2 - z_{k-1}^2),\ \tfrac13(z_k^3 - z_{k-1}^3)\right)

    Parameters
    ----------
    E1, E2, nu12, G12 : float
        Ply properties, as in :func:`lamina_stiffness`.
    angles : array_like
        Fibre angle of each ply [rad], listed from the bottom (most negative
        z) to the top.
    ply_thickness : float or array_like
        Thickness of every ply, or of each [m], positive.

    Returns
    -------
    Laminate
        The stacking sequence, the A, B and D matrices and the midplane
        and ply-stress solvers.

    Notes
    -----
    Kirchhoff plate kinematics (straight normals) and perfect bonding are
    assumed; interlaminar stresses near free edges are not modelled. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    E. Reissner and Y. Stavsky, "Bending and stretching of certain types of
    heterogeneous aeolotropic elastic plates," *J. Appl. Mech.* 28 (1961)
    402-408.

    Examples
    --------
    >>> lam = classical_laminate(140e9, 10e9, 0.3, 5e9, [0.0, np.pi / 2, np.pi / 2, 0.0], 0.125e-3)
    >>> bool(np.allclose(lam.B, 0))
    True
    """
    Q = lamina_stiffness(E1, E2, nu12, G12)
    angles = _vector("angles", angles)
    t = np.broadcast_to(np.asarray(ply_thickness, float), angles.shape)
    _positive(ply_thickness=t)
    z = np.concatenate([[0.0], np.cumsum(t)]) - t.sum() / 2
    Qbar = np.array([transformed_stiffness(Q, a) for a in angles])
    A = np.einsum("kij,k->ij", Qbar, np.diff(z))
    B = np.einsum("kij,k->ij", Qbar, np.diff(z**2) / 2)
    D = np.einsum("kij,k->ij", Qbar, np.diff(z**3) / 3)
    return Laminate(angles=angles, z=z, Q=Q, Qbar=Qbar, A=A, B=B, D=D)
