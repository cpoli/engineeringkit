r"""Machine elements: power screws, journal bearings and rolling bearings.

Coulomb's friction law turns a screw thread into an inclined plane wrapped
round a cylinder. Petroff treated a lightly loaded journal bearing as viscous
shear across a thin concentric film; Reynolds derived the pressure in an
eccentric film and Sommerfeld solved it for a long bearing. Lundberg and
Palmgren gave the rolling-bearing life :math:`L_{10} = (C/P)^p`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.mechanical_engineering.utils.roots import bracketed_roots

__all__ = [
    "screw_torque_raise",
    "petroff_friction_torque",
    "sommerfeld_number",
    "long_bearing_eccentricity",
    "bearing_l10_life",
]


def screw_torque_raise(load, mean_diameter, lead, mu):
    r"""Torque to raise a load with a square-thread power screw.

    .. math:: T = \frac{W d_m}{2}\tan(\phi + \lambda), \quad
              \tan\lambda = \frac{L}{\pi d_m}, \quad \tan\phi = \mu

    Parameters
    ----------
    load : float or array_like
        Axial load :math:`W` [N], nonnegative.
    mean_diameter : float or array_like
        Mean thread diameter :math:`d_m` [m], positive.
    lead : float or array_like
        Axial advance per revolution :math:`L` [m], positive.
    mu : float or array_like
        Thread friction coefficient [-], nonnegative.

    Returns
    -------
    float or ndarray
        Torque [N m], excluding collar friction.

    Raises
    ------
    ValueError
        If :math:`\phi + \lambda \ge \pi/2`, where no torque can raise the load.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    C. A. Coulomb, "Théorie des machines simples," *Mém. Math. Phys. Acad.
    R. Sci.* 10 (1785) 161-332.

    Examples
    --------
    >>> round(float(screw_torque_raise(1000.0, 0.02, 0.005, 0.0) * 2 * np.pi), 6)
    5.0
    """
    _finite(load=load, mean_diameter=mean_diameter, lead=lead, mu=mu)
    _nonnegative(load=load, mu=mu)
    _positive(mean_diameter=mean_diameter, lead=lead)
    lam = np.arctan(lead / (np.pi * mean_diameter))
    phi = np.arctan(mu)
    if np.any(phi + lam >= np.pi / 2):
        raise ValueError("friction and lead angles must sum to less than pi/2")
    return np.asarray(load) * mean_diameter / 2 * np.tan(phi + lam)


def petroff_friction_torque(viscosity, radius, length, clearance, omega):
    r"""Petroff's viscous friction torque in a concentric journal bearing.

    .. math:: T = \frac{2\pi \mu \omega r^3 L}{c}

    Parameters
    ----------
    viscosity : float or array_like
        Dynamic viscosity of the lubricant :math:`\mu` [Pa s], positive.
    radius : float or array_like
        Journal radius :math:`r` [m], positive.
    length : float or array_like
        Bearing length :math:`L` [m], positive.
    clearance : float or array_like
        Radial clearance :math:`c` [m], positive and much less than ``radius``.
    omega : float or array_like
        Journal speed [rad/s].

    Returns
    -------
    float or ndarray
        Friction torque [N m], with the sign of ``omega``.

    Notes
    -----
    The film is a linear (Couette) shear flow across a gap :math:`c \ll r`
    with the journal centred, which holds for light loads and high speeds.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    N. P. Petroff, "Friction in machines and the effect of the lubricant,"
    *Inzhenernyi Zhurnal* (St Petersburg), nos. 1-4 (1883), in Russian.

    Examples
    --------
    >>> round(float(petroff_friction_torque(0.05, 0.02, 0.04, 2e-5, 100.0)), 6)
    0.502655
    """
    _finite(omega=omega)
    _positive(viscosity=viscosity, radius=radius, length=length, clearance=clearance)
    return 2 * np.pi * viscosity * omega * radius**3 * length / clearance


def sommerfeld_number(viscosity, omega, load, radius, length, clearance):
    r"""Bearing characteristic number :math:`S = (r/c)^2 \mu N / P`.

    Parameters
    ----------
    viscosity : float or array_like
        Dynamic viscosity :math:`\mu` [Pa s], positive.
    omega : float or array_like
        Journal speed [rad/s], positive; :math:`N = \omega/2\pi` [rev/s].
    load : float or array_like
        Radial load :math:`W` [N], positive.
    radius, length, clearance : float or array_like
        Journal radius, bearing length and radial clearance [m], positive.

    Returns
    -------
    float or ndarray
        Sommerfeld number [-], with projected pressure :math:`P = W/(2rL)`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    A. Sommerfeld, "Zur hydrodynamischen Theorie der Schmiermittelreibung,"
    *Z. Math. Phys.* 50 (1904) 97-155.

    Examples
    --------
    >>> round(float(sommerfeld_number(0.05, 2 * np.pi * 30, 2000.0, 0.02, 0.04, 2e-5)), 4)
    1.2
    """
    _positive(viscosity=viscosity, omega=omega, load=load, radius=radius, length=length, clearance=clearance)
    pressure = load / (2 * radius * length)
    return (radius / clearance) ** 2 * viscosity * (omega / (2 * np.pi)) / pressure


def _sommerfeld_full(eps):
    return (2 + eps * eps) * np.sqrt(1 - eps * eps) / (12 * np.pi**2 * eps)


def long_bearing_eccentricity(S):
    r"""Eccentricity ratio of an infinitely long journal bearing (full Sommerfeld solution).

    Solves

    .. math:: S = \frac{(2 + \varepsilon^2)\sqrt{1 - \varepsilon^2}}{12\pi^2\varepsilon}

    for :math:`\varepsilon = e/c`, the journal's offset as a fraction of the
    radial clearance.

    Parameters
    ----------
    S : float or array_like
        Sommerfeld number [-], positive (see :func:`sommerfeld_number`).

    Returns
    -------
    float or ndarray
        Eccentricity ratio [-] in (0, 1); it tends to zero (Petroff's
        centred journal) as :math:`S \to \infty`.

    Notes
    -----
    The full Sommerfeld solution keeps the negative film pressures that a
    real oil film cavitates away, so it puts the journal at 90 degrees to
    the load. The short-bearing and finite-bearing solutions are not
    implemented. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    O. Reynolds, "On the theory of lubrication and its application to Mr.
    Beauchamp Tower's experiments," *Phil. Trans. R. Soc. Lond.* 177 (1886)
    157-234; Sommerfeld (1904), as for :func:`sommerfeld_number`.

    Examples
    --------
    >>> eps = long_bearing_eccentricity(0.1)
    >>> round(float(eps), 4)
    0.1688
    """
    _positive(S=S)
    grid = np.concatenate((np.geomspace(1e-12, 0.5, 400), 1 - np.geomspace(0.5, 1e-15, 400)[1:]))
    s = np.asarray(S, dtype=float)
    eps = np.array([bracketed_roots(lambda e, target=value: _sommerfeld_full(e) - target, grid)[0] for value in s.ravel()])
    return eps.reshape(s.shape) if s.ndim else float(eps[0])


def bearing_l10_life(C, P, p=3):
    r"""Lundberg-Palmgren rating life :math:`L_{10} = (C/P)^p` of a rolling bearing.

    Parameters
    ----------
    C : float or array_like
        Basic dynamic load rating [N], positive.
    P : float or array_like
        Equivalent dynamic load [N], positive.
    p : float, optional
        Life exponent [-], positive: 3 for ball bearings, 10/3 for roller
        bearings.

    Returns
    -------
    float or ndarray
        Revolutions that 90% of a large group of identical bearings complete
        before the first sign of fatigue.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    G. Lundberg and A. Palmgren, "Dynamic capacity of rolling bearings,"
    *Acta Polytech. Mech. Eng. Ser.* 1, no. 3 (1947).

    Examples
    --------
    >>> bearing_l10_life(2.0, 1.0)
    8000000.0
    """
    _finite(C=C, P=P, p=p)
    _positive(C=C, P=P, p=p)
    return (C / P) ** p * 1e6
