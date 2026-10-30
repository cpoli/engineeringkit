r"""Column buckling and design: Euler, Engesser, the secant formula and Johnson's parabola.

An ideal pin-ended column of length :math:`L` buckles at
:math:`P_{cr} = \pi^2EI/(KL)^2`, where the effective-length factor :math:`K`
accounts for the end conditions. When the buckling stress exceeds the
proportional limit, Engesser replaced :math:`E` by the tangent modulus
:math:`E_t = d\sigma/d\varepsilon` at the buckling stress. Real columns are
loaded off-centre; Young's secant formula gives the peak stress of an
eccentrically loaded column, and Johnson's parabola replaces Euler's
hyperbola by a design curve where it would exceed half the yield stress.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "euler_buckling_load",
    "radius_of_gyration",
    "slenderness_ratio",
    "tangent_modulus_load",
    "effective_length_factor",
    "secant_formula_stress",
    "secant_formula_load",
    "johnson_critical_stress",
]

_END_CONDITIONS = ("pinned-pinned", "fixed-pinned", "fixed-fixed", "fixed-free", "fixed-guided")


def euler_buckling_load(E, I, L, K=1.0):
    r"""Elastic critical load of an ideal column, :math:`P_{cr} = \pi^2EI/(KL)^2`.

    Parameters
    ----------
    E : float or array_like
        Young's modulus [Pa], positive.
    I : float or array_like
        Least second moment of area [m^4], positive.
    L : float or array_like
        Column length [m], positive.
    K : float or array_like, optional
        Effective-length factor [-], positive: 1 pinned-pinned, 0.5
        fixed-fixed, about 0.7 fixed-pinned, 2 fixed-free.

    Returns
    -------
    float or ndarray
        Critical load [N].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    L. Euler, *Methodus inveniendi lineas curvas maximi minimive proprietate
    gaudentes*, Additamentum I, "De curvis elasticis" (Lausanne, 1744).
    """
    _finite(E=E, I=I, L=L, K=K)
    _positive(E=E, I=I, L=L, K=K)
    return np.pi**2 * np.asarray(E) * I / (np.asarray(K) * L) ** 2


def radius_of_gyration(I, A):
    r"""Radius of gyration, :math:`r = \sqrt{I/A}`.

    Parameters
    ----------
    I : float or array_like
        Second moment of area [m^4], positive.
    A : float or array_like
        Area [m^2], positive.

    Returns
    -------
    float or ndarray
        Radius of gyration [m].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(I=I, A=A)
    _positive(I=I, A=A)
    return np.sqrt(np.asarray(I) / A)


def slenderness_ratio(L, K, r):
    r"""Effective slenderness ratio, :math:`KL/r`.

    Parameters
    ----------
    L : float or array_like
        Column length [m], positive.
    K : float or array_like
        Effective-length factor [-], positive.
    r : float or array_like
        Radius of gyration [m], positive.

    Returns
    -------
    float or ndarray
        Slenderness [-]; the Euler stress is :math:`\pi^2E/(KL/r)^2`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(L=L, K=K, r=r)
    _positive(L=L, K=K, r=r)
    return np.asarray(K) * L / r


def tangent_modulus_load(Et, I, L, K=1.0, *, area=None):
    r"""Engesser's tangent-modulus buckling load, :math:`P_t = \pi^2E_tI/(KL)^2`.

    If ``Et`` is a number, the formula is evaluated directly. If it is a
    function of stress, the buckling stress must satisfy

    .. math::

        \sigma = \frac{\pi^2E_t(\sigma)}{(KL/r)^2},\qquad r^2 = I/A,

    which is solved for :math:`\sigma` by Brent's method between zero and the
    Euler stress computed with the initial modulus :math:`E_t(0)`.

    Parameters
    ----------
    Et : float, array_like or callable
        Tangent modulus [Pa], positive, or a function returning the tangent
        modulus [Pa] at a compressive stress [Pa]. A function must be
        positive at zero stress, nonnegative above it, and must not
        increase with stress.
    I : float or array_like
        Least second moment of area [m^4], positive.
    L : float or array_like
        Column length [m], positive.
    K : float or array_like, optional
        Effective-length factor [-], positive.
    area : float, optional
        Cross-sectional area [m^2], positive. Required when ``Et`` is a
        function; ``I``, ``L`` and ``K`` must then be scalars.

    Returns
    -------
    float or ndarray
        Critical load [N]. It equals :func:`euler_buckling_load` wherever the
        column buckles while still elastic.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    F. Engesser, "Ueber die Knickfestigkeit gerader Stäbe," *Z. Architekten-
    und Ingenieur-Vereins zu Hannover* 35 (1889) 455-462; F. R. Shanley,
    "Inelastic column theory," *J. Aeronaut. Sci.* 14 (1947) 261-268.
    """
    if not callable(Et):
        return euler_buckling_load(Et, I, L, K)
    if area is None:
        raise ValueError("area is required when Et is a function of stress")
    _positive(I=I, L=L, K=K, area=area)
    if np.ndim(I) or np.ndim(L) or np.ndim(K) or np.ndim(area):
        raise ValueError("I, L, K and area must be scalars when Et is a function of stress")
    ratio = np.pi**2 * I / (area * (K * L) ** 2)

    def modulus(sigma):
        value = float(Et(sigma))
        _nonnegative(Et=value)
        return value

    initial = modulus(0.0)
    _positive(Et=initial)
    upper = ratio * initial
    if upper - ratio * modulus(upper) < 0:
        raise ValueError("Et must not increase with stress")
    sigma = brentq(lambda s: s - ratio * modulus(s), 0.0, upper, xtol=1e-12 * upper, rtol=1e-14)
    return float(sigma) * area


def effective_length_factor(ends):
    r"""Theoretical effective-length factor :math:`K` of an ideal column.

    :math:`KL` is the length of the pin-ended column with the same Euler
    load, the distance between inflection points of the buckled shape. For
    a column fixed at one end and pinned at the other, the buckling
    condition is :math:`\tan kL = kL`, whose first root
    :math:`kL = 4.4934` gives :math:`K = \pi/4.4934 = 0.699`.

    Parameters
    ----------
    ends : {"pinned-pinned", "fixed-pinned", "fixed-fixed", "fixed-free", "fixed-guided"}
        End conditions; ``"fixed-guided"`` has both ends restrained against
        rotation but one free to sway.

    Returns
    -------
    float
        :math:`K` [-]: 1, 0.699, 0.5, 2 and 1 respectively.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Design codes recommend larger values (0.65 for
    fixed-fixed, 0.80 for fixed-pinned, 2.1 for fixed-free) because real
    connections are never perfectly rigid.
    """
    if ends not in _END_CONDITIONS:
        raise ValueError(f"ends must be one of {_END_CONDITIONS}")
    if ends == "fixed-pinned":
        root = brentq(lambda kL: np.sin(kL) - kL * np.cos(kL), np.pi, 1.5 * np.pi, xtol=1e-14)
        return float(np.pi / root)
    return {"pinned-pinned": 1.0, "fixed-fixed": 0.5, "fixed-free": 2.0, "fixed-guided": 1.0}[ends]


def secant_formula_stress(P, A, e, c, r, L, E, K=1.0):
    r"""Peak stress in an eccentrically loaded column, the secant formula.

    .. math::

        \sigma_\max = \frac{P}{A}\left[1 + \frac{ec}{r^2}
        \sec\left(\frac{KL}{2r}\sqrt{\frac{P}{EA}}\right)\right].

    The stress is not proportional to the load: it grows without bound as
    :math:`P` approaches the Euler load.

    Parameters
    ----------
    P : float or array_like
        Compressive load [N], positive and below the Euler load.
    A : float or array_like
        Area [m^2], positive.
    e : float or array_like
        Eccentricity of the load [m], nonnegative.
    c : float or array_like
        Distance from the centroid to the extreme compressive fibre [m], positive.
    r : float or array_like
        Radius of gyration in the plane of bending [m], positive.
    L : float or array_like
        Length [m], positive.
    E : float or array_like
        Young's modulus [Pa], positive.
    K : float or array_like, optional
        Effective-length factor [-], positive.

    Returns
    -------
    float or ndarray
        Largest compressive stress [Pa], as a positive number.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Linear elastic: valid while the stress stays
    below the proportional limit.

    References
    ----------
    T. Young, *A Course of Lectures on Natural Philosophy and the Mechanical
    Arts* (London: Johnson, 1807), vol. 2.
    """
    _finite(P=P, A=A, e=e, c=c, r=r, L=L, E=E, K=K)
    _positive(P=P, A=A, c=c, r=r, L=L, E=E, K=K)
    _nonnegative(e=e)
    angle = np.asarray(K) * L / (2 * np.asarray(r)) * np.sqrt(np.asarray(P) / (np.asarray(E) * A))
    if np.any(angle >= np.pi / 2):
        raise ValueError("P must be below the Euler load")
    return np.asarray(P) / A * (1 + np.asarray(e) * c / np.asarray(r) ** 2 / np.cos(angle))


def secant_formula_load(sigma_y, A, e, c, r, L, E, K=1.0):
    r"""Load at which the secant-formula stress reaches yield.

    Solves :func:`secant_formula_stress` :math:`= \sigma_y` for :math:`P` by
    Brent's method. With :math:`e = 0` it is the smaller of the squash load
    :math:`\sigma_yA` and the Euler load.

    Parameters
    ----------
    sigma_y : float
        Yield stress [Pa], positive.
    A : float
        Area [m^2], positive.
    e : float
        Eccentricity [m], nonnegative.
    c : float
        Distance from the centroid to the extreme compressive fibre [m], positive.
    r : float
        Radius of gyration [m], positive.
    L : float
        Length [m], positive.
    E : float
        Young's modulus [Pa], positive.
    K : float, optional
        Effective-length factor [-], positive.

    Returns
    -------
    float
        Load at first yield [N].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Scalars only.
    """
    _positive(sigma_y=sigma_y, A=A, c=c, r=r, L=L, E=E, K=K)
    _nonnegative(e=e)
    if any(np.ndim(v) for v in (sigma_y, A, e, c, r, L, E, K)):
        raise ValueError("secant_formula_load takes scalars")
    euler = float(euler_buckling_load(E, A * r**2, L, K))
    squash = sigma_y * A
    if e == 0:
        return min(squash, euler)
    ratio = e * c / r**2

    def excess(P):
        cosine = np.cos(K * L / (2 * r) * np.sqrt(P / (E * A)))
        return P / A * (cosine + ratio) - sigma_y * cosine

    upper = min(squash, euler)
    return float(brentq(excess, 0.0, upper, xtol=1e-12 * upper, rtol=1e-14))


def johnson_critical_stress(sigma_y, E, slenderness):
    r"""Critical stress on Johnson's parabola, joined to Euler's hyperbola.

    .. math::

        \sigma_{cr} = \sigma_y - \frac{\sigma_y^2}{4\pi^2E}\left(\frac{KL}{r}\right)^2
        \quad\text{for}\quad \frac{KL}{r} < C_c = \sqrt{\frac{2\pi^2E}{\sigma_y}},

    and :math:`\pi^2E/(KL/r)^2` beyond. The parabola starts at the yield
    stress for a stub column and meets Euler's curve tangentially at
    :math:`\sigma_y/2`.

    Parameters
    ----------
    sigma_y : float or array_like
        Yield stress [Pa], positive.
    E : float or array_like
        Young's modulus [Pa], positive.
    slenderness : float or array_like
        Effective slenderness ratio :math:`KL/r` [-], nonnegative.

    Returns
    -------
    float or ndarray
        Critical stress [Pa].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. The curve includes no safety factor.

    References
    ----------
    J. B. Johnson, C. W. Bryan and F. E. Turneaure, *The Theory and Practice
    of Modern Framed Structures* (New York: Wiley, 1893).
    """
    _positive(sigma_y=sigma_y, E=E)
    _nonnegative(slenderness=slenderness)
    sy, s = np.asarray(sigma_y, float), np.asarray(slenderness, float)
    transition = np.sqrt(2 * np.pi**2 * np.asarray(E) / sy)
    euler = np.pi**2 * np.asarray(E) / np.where(s > 0, s, 1.0) ** 2
    return np.where(s < transition, sy - sy**2 * s**2 / (4 * np.pi**2 * np.asarray(E)), euler)
