r"""Reliability of systems built from independent components.

A series system works only if every component works; a parallel system
needs one; a :math:`k`-out-of-:math:`n` system needs any :math:`k`. Von
Neumann showed that redundancy can make a reliable machine from unreliable
parts, and Birnbaum measured how much each part matters. Every formula here
assumes independent component failures.
"""

from __future__ import annotations

import math

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval

__all__ = ["series_reliability", "parallel_reliability", "k_of_n_reliability", "birnbaum_importance_series"]


def series_reliability(*r):
    r"""Reliability of independent components in series, :math:`R_s = \prod_i R_i`.

    Parameters
    ----------
    *r : float
        Component reliabilities [-], each in [0, 1].

    Returns
    -------
    float
        System reliability [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(series_reliability(0.9, 0.8), 2)
    0.72
    """
    _interval("r", r, 0, 1)
    return float(np.prod(r))


def parallel_reliability(*r):
    r"""Reliability of independent components in active parallel, :math:`1 - \prod_i (1 - R_i)`.

    Parameters
    ----------
    *r : float
        Component reliabilities [-], each in [0, 1].

    Returns
    -------
    float
        System reliability [-].

    Notes
    -----
    Assumes perfect switching and no common-cause failures. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    J. von Neumann, "Probabilistic logics and the synthesis of reliable
    organisms from unreliable components," in *Automata Studies* (Princeton,
    1956) 43-98.

    Examples
    --------
    >>> round(parallel_reliability(0.9, 0.8), 2)
    0.98
    """
    _interval("r", r, 0, 1)
    return 1 - float(np.prod([1 - x for x in r]))


def k_of_n_reliability(n, k, r):
    r"""Reliability of a :math:`k`-out-of-:math:`n` system of identical independent components.

    :math:`R = \sum_{i=k}^{n} \binom{n}{i} r^i (1-r)^{n-i}`.

    Parameters
    ----------
    n : int
        Number of components, at least one.
    k : int
        Number that must work, from zero to ``n``.
    r : float
        Component reliability [-], in [0, 1].

    Returns
    -------
    float
        System reliability [-].

    Notes
    -----
    ``k = n`` is a series system and ``k = 1`` a parallel one. Inputs are raw
    SI numbers with no dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    Z. W. Birnbaum, J. D. Esary and S. C. Saunders, "Multi-component systems
    and structures and their reliability," *Technometrics* 3 (1961) 55-77.

    Examples
    --------
    >>> round(k_of_n_reliability(3, 2, 0.9), 3)
    0.972
    """
    _finite(n=n, k=k, r=r)
    _integer("n", n, 1)
    _integer("k", k, 0)
    _interval("r", r, 0, 1)
    if k > n:
        raise ValueError("k cannot exceed n")
    return sum(math.comb(n, i) * r**i * (1 - r) ** (n - i) for i in range(k, n + 1))


def birnbaum_importance_series(reliabilities, i):
    r"""Birnbaum importance :math:`\partial R_s/\partial R_i = \prod_{j \ne i} R_j` of a series component.

    Parameters
    ----------
    reliabilities : array_like
        Component reliabilities [-], each in [0, 1].
    i : int
        Index of the component, from zero.

    Returns
    -------
    float
        Sensitivity of system reliability to component ``i`` [-].

    Notes
    -----
    In a series system the *least* reliable component has the largest
    importance. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    Z. W. Birnbaum, "On the importance of different components in a
    multicomponent system," in *Multivariate Analysis II* (Academic Press,
    1969) 581-592.

    Examples
    --------
    >>> round(birnbaum_importance_series([0.7, 0.8, 0.9], 0), 2)
    0.72
    """
    _finite(reliabilities=reliabilities, i=i)
    _interval("reliabilities", reliabilities, 0, 1)
    _integer("i", i, 0)
    if i >= len(reliabilities):
        raise ValueError("component index out of range")
    return float(np.prod(np.delete(np.asarray(reliabilities), i)))
