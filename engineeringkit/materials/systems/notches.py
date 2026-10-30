r"""Notches: stress concentration, notch sensitivity and Neuber's rule.

A hole or a fillet raises the elastic stress at its root by the factor
:math:`K_t`; Inglis (1913) found :math:`1 + 2a/b` for an elliptical hole. In
fatigue a sharp notch hurts less than :math:`K_t` suggests, because the
stress falls off steeply within a small material length; the notch
sensitivity :math:`q` gives the fatigue factor :math:`K_f = 1 + q(K_t - 1)`.
Once the root yields, Neuber (1961) showed that the product of the local
stress and strain concentrations stays :math:`K_t^2`, so the elastic
solution still fixes the local stress-strain state.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit.materials.systems.plasticity import ramberg_osgood_strain

__all__ = [
    "elliptical_hole_concentration",
    "peterson_notch_sensitivity",
    "neuber_notch_sensitivity",
    "fatigue_notch_factor",
    "neuber_notch_stress",
]


def elliptical_hole_concentration(a, b):
    r"""Inglis's stress concentration at an elliptical hole in a wide plate, :math:`K_t = 1 + 2a/b`.

    Parameters
    ----------
    a : float or array_like
        Semi-axis across the load [m], positive.
    b : float or array_like
        Semi-axis along the load [m], positive.

    Returns
    -------
    float or ndarray
        :math:`K_t` [-]: 3 for a circular hole, unbounded as the ellipse
        closes into a crack.

    Notes
    -----
    Equivalently :math:`K_t = 1 + 2\sqrt{a/\rho}` with the root radius
    :math:`\rho = b^2/a`. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/materials`.

    References
    ----------
    C. E. Inglis, "Stresses in a plate due to the presence of cracks and
    sharp corners," *Trans. Inst. Naval Archit.* 55 (1913) 219-241.

    Examples
    --------
    >>> float(elliptical_hole_concentration(1.0, 1.0))
    3.0
    """
    _positive(a=a, b=b)
    return 1 + 2 * np.asarray(a) / b


def peterson_notch_sensitivity(radius, alpha):
    r"""Peterson's notch sensitivity :math:`q = 1/(1 + \alpha/r)`.

    Parameters
    ----------
    radius : float or array_like
        Notch root radius :math:`r` [m], positive.
    alpha : float or array_like
        Material characteristic length [m], positive (about 0.25 mm for
        normalised steels, falling with strength).

    Returns
    -------
    float or ndarray
        Notch sensitivity [-], between 0 and 1.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    R. E. Peterson, "Notch-sensitivity," in G. Sines and J. L. Waisman
    (eds.), *Metal Fatigue* (McGraw-Hill, New York, 1959), 293-306.
    """
    _positive(radius=radius, alpha=alpha)
    return 1 / (1 + np.asarray(alpha) / radius)


def neuber_notch_sensitivity(radius, rho):
    r"""Neuber's notch sensitivity :math:`q = 1/(1 + \sqrt{\rho/r})`.

    Parameters
    ----------
    radius : float or array_like
        Notch root radius :math:`r` [m], positive.
    rho : float or array_like
        Neuber's material length :math:`\rho` [m], positive.

    Returns
    -------
    float or ndarray
        Notch sensitivity [-], between 0 and 1.

    Notes
    -----
    Shigley tabulates :math:`\sqrt{\rho}` against ultimate strength.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _positive(radius=radius, rho=rho)
    return 1 / (1 + np.sqrt(np.asarray(rho) / radius))


def fatigue_notch_factor(Kt, q):
    r"""Fatigue stress concentration :math:`K_f = 1 + q(K_t - 1)`.

    Parameters
    ----------
    Kt : float or array_like
        Elastic stress concentration factor [-], at least 1.
    q : float or array_like
        Notch sensitivity [-], in [0, 1].

    Returns
    -------
    float or ndarray
        :math:`K_f` [-], between 1 and :math:`K_t`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> float(fatigue_notch_factor(3.0, 0.5))
    2.0
    """
    _interval("Kt", Kt, 1, np.inf)
    _interval("q", q, 0, 1)
    return 1 + np.asarray(q) * (np.asarray(Kt) - 1)


def neuber_notch_stress(nominal_stress, Kt, E, sigma_y, n, offset=0.002):
    r"""Local notch-root stress and strain by Neuber's rule on a Ramberg-Osgood curve.

    Solves :math:`\sigma\varepsilon = (K_t S)^2/E` together with
    :math:`\varepsilon = \sigma/E + \varepsilon_0(\sigma/\sigma_y)^n`.

    Parameters
    ----------
    nominal_stress : float or array_like
        Nominal (net-section) stress :math:`S` [Pa], elastic; its sign is
        kept.
    Kt : float or array_like
        Elastic stress concentration factor [-], at least 1 (use
        :math:`K_f` for fatigue).
    E : float or array_like
        Young's modulus [Pa], positive.
    sigma_y : float or array_like
        Offset yield strength of the (cyclic) curve [Pa], positive.
    n : float or array_like
        Ramberg-Osgood hardening exponent [-], at least 1.
    offset : float, optional
        Plastic strain at ``sigma_y`` [-], positive.

    Returns
    -------
    tuple of float or ndarray
        ``(stress, strain)`` at the notch root [Pa], [-]. The stress is at
        most :math:`K_t S` and the strain at least :math:`K_t S/E`.

    Notes
    -----
    Newton's method from the elastic stress :math:`K_t S` converges
    monotonically, since :math:`\sigma\varepsilon(\sigma)` is increasing and
    convex. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    H. Neuber, "Theory of stress concentration for shear-strained
    prismatical bodies with arbitrary nonlinear stress-strain law," *J.
    Appl. Mech.* 28 (1961) 544-550.

    Examples
    --------
    >>> s, e = neuber_notch_stress(100e6, 2.0, 200e9, 1e9, 10)
    >>> round(float(s) / 1e6, 3)
    200.0
    """
    _finite(nominal_stress=nominal_stress)
    _interval("Kt", Kt, 1, np.inf)
    _positive(E=E, sigma_y=sigma_y, offset=offset)
    _interval("n", n, 1, np.inf)
    S = np.asarray(nominal_stress, float)
    elastic = np.asarray(Kt) * np.abs(S)
    target = elastic * elastic / E
    # sigma = 0 only where S = 0, which is already the root; dividing by 1 there keeps the step zero.
    sigma = elastic * np.ones(np.broadcast(elastic, E, sigma_y, n, offset).shape)
    for _ in range(100):
        eps = ramberg_osgood_strain(sigma, E, sigma_y, n, offset)
        plastic = eps - sigma / E
        slope = 2 * sigma / E + plastic * (n + 1)
        step = (sigma * eps - target) / np.where(slope > 0, slope, 1)
        sigma = sigma - step
        if np.all(np.abs(step) <= 1e-12 * np.asarray(sigma_y)):
            break
    strain = ramberg_osgood_strain(sigma, E, sigma_y, n, offset)
    sign = np.where(S < 0, -1.0, 1.0)
    return sign * sigma, sign * strain
