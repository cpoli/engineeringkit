r"""Plasticity: dislocations, grain-size strengthening and the Ramberg-Osgood curve.

Frenkel estimated that shearing a perfect crystal needs about :math:`G/2\pi`,
thousands of times what real metals need. Taylor, Orowan and Polanyi
explained the gap with the dislocation, a line defect that moves one row of
atoms at a time against the small Peierls stress. Hall and Petch found that
grain boundaries block dislocations, so yield strength grows as
:math:`d^{-1/2}`. Ramberg and Osgood described the smooth curve from elastic
to plastic by three parameters.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["ideal_shear_strength", "peierls_stress", "hall_petch", "ramberg_osgood_strain"]


def ideal_shear_strength(G, b, a):
    r"""Frenkel's theoretical shear strength of a perfect crystal.

    .. math:: \tau_{max} = \frac{G\,b}{2\pi a}

    Parameters
    ----------
    G : float or array_like
        Shear modulus [Pa], positive.
    b : float or array_like
        Atomic spacing in the slip direction [m], positive.
    a : float or array_like
        Spacing of the slip planes [m], positive.

    Returns
    -------
    float or ndarray
        Stress to slide one plane rigidly over the next [Pa].

    Notes
    -----
    A sinusoidal resistance whose initial slope is :math:`G`. Inputs are raw
    SI numbers with no dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    J. Frenkel, "Zur Theorie der Elastizitätsgrenze und der Festigkeit
    kristallinischer Körper," *Z. Phys.* 37 (1926) 572-609.

    Examples
    --------
    >>> round(float(ideal_shear_strength(2 * np.pi, 1.0, 1.0)), 12)
    1.0
    """
    _positive(G=G, b=b, a=a)
    return np.asarray(G) * b / (2 * np.pi * np.asarray(a))


def peierls_stress(G, nu, b, a):
    r"""Peierls-Nabarro stress to move an edge dislocation through a perfect lattice.

    .. math:: \tau_P = \frac{2G}{1 - \nu}\exp\!\left(-\frac{2\pi a}{(1 - \nu)\,b}\right)

    Parameters
    ----------
    G : float or array_like
        Shear modulus [Pa], positive.
    nu : float or array_like
        Poisson's ratio [-], in (-1, 0.5).
    b : float or array_like
        Burgers vector [m], positive.
    a : float or array_like
        Slip-plane spacing [m], positive.

    Returns
    -------
    float or ndarray
        Lattice friction stress [Pa].

    Notes
    -----
    The exponential makes the estimate a guide to trends (wide, close-packed
    planes slip easily) rather than an accurate value. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    R. Peierls, "The size of a dislocation," *Proc. Phys. Soc.* 52 (1940)
    34-37; F. R. N. Nabarro, "Dislocations in a simple cubic lattice,"
    *Proc. Phys. Soc.* 59 (1947) 256-272.
    """
    _positive(G=G, b=b, a=a)
    _interval("nu", nu, -1, 0.5, closed=False)
    nu = np.asarray(nu)
    return 2 * np.asarray(G) / (1 - nu) * np.exp(-2 * np.pi * np.asarray(a) / ((1 - nu) * b))


def hall_petch(sigma0, k, d):
    r"""Hall-Petch yield strength :math:`\sigma_y = \sigma_0 + k\,d^{-1/2}`.

    Parameters
    ----------
    sigma0 : float or array_like
        Friction stress of a single crystal [Pa], nonnegative.
    k : float or array_like
        Strengthening coefficient [Pa sqrt(m)], nonnegative.
    d : float or array_like
        Mean grain diameter [m], positive.

    Returns
    -------
    float or ndarray
        Yield strength [Pa].

    Notes
    -----
    Holds down to grain sizes of tens of nanometres, below which grain
    boundary sliding softens the metal (the inverse Hall-Petch effect).
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    E. O. Hall, "The deformation and ageing of mild steel: III Discussion of
    results," *Proc. Phys. Soc. B* 64 (1951) 747-753; N. J. Petch, "The
    cleavage strength of polycrystals," *J. Iron Steel Inst.* 174 (1953)
    25-28.
    """
    _nonnegative(sigma0=sigma0, k=k)
    _positive(d=d)
    return sigma0 + k / np.sqrt(np.asarray(d))


def ramberg_osgood_strain(stress, E, yield_strength, n, offset=0.002):
    r"""Total strain on a Ramberg-Osgood stress-strain curve.

    .. math:: \varepsilon = \frac{\sigma}{E} + \varepsilon_o\left(\frac{\sigma}{\sigma_y}\right)^n

    Parameters
    ----------
    stress : float or array_like
        Stress [Pa]; the curve is odd, so compression mirrors tension.
    E : float or array_like
        Young's modulus [Pa], positive.
    yield_strength : float or array_like
        Offset yield strength :math:`\sigma_y` [Pa], positive.
    n : float or array_like
        Hardening exponent [-], at least 1 (large ``n`` is nearly
        elastic-perfectly plastic).
    offset : float or array_like
        Plastic strain :math:`\varepsilon_o` at :math:`\sigma_y` [-],
        positive; 0.002 is the 0.2% proof stress.

    Returns
    -------
    float or ndarray
        Strain [-], elastic plus plastic.

    Notes
    -----
    Ramberg and Osgood wrote :math:`\varepsilon = \sigma/E + K(\sigma/E)^n`;
    the offset form is the same curve with :math:`K = \varepsilon_o (E/\sigma_y)^n`.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    W. Ramberg and W. R. Osgood, "Description of stress-strain curves by
    three parameters," NACA Technical Note 902 (1943).

    Examples
    --------
    >>> round(float(ramberg_osgood_strain(400e6, 200e9, 400e6, 10)), 6)
    0.004
    """
    _finite(stress=stress)
    _positive(E=E, yield_strength=yield_strength, offset=offset)
    _interval("n", n, 1, np.inf, closed="left")
    s = np.asarray(stress, float)
    return s / E + offset * np.sign(s) * (np.abs(s) / yield_strength) ** n
