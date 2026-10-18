r"""Column buckling: Euler's elastic load and Engesser's tangent-modulus load.

An ideal pin-ended column of length :math:`L` buckles at
:math:`P_{cr} = \pi^2EI/(KL)^2`, where the effective-length factor :math:`K`
accounts for the end conditions. When the buckling stress exceeds the
proportional limit, Engesser replaced :math:`E` by the tangent modulus
:math:`E_t = d\sigma/d\varepsilon` at the buckling stress.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["euler_buckling_load", "radius_of_gyration", "slenderness_ratio", "tangent_modulus_load"]


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
