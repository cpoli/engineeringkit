r"""Reliability of systems built from independent components.

A series system works only if every component works; a parallel system
needs one; a :math:`k`-out-of-:math:`n` system needs any :math:`k`. Von
Neumann showed that redundancy can make a reliable machine from unreliable
parts, and Birnbaum measured how much each part matters. A general block
diagram, such as a bridge, is a network whose minimal path and cut sets
give its exact reliability and Esary and Proschan's bounds. A cold standby
unit is switched in only when the running one fails. Every formula here
assumes independent component failures.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "series_reliability",
    "parallel_reliability",
    "k_of_n_reliability",
    "birnbaum_importance_series",
    "minimal_path_sets",
    "minimal_cut_sets",
    "network_reliability",
    "esary_proschan_bounds",
    "standby_reliability",
    "standby_mttf",
]

_MAX_COMPONENTS = 20


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
    parts = np.broadcast_arrays(*(np.asarray(x) for x in r))
    _interval("r", np.array(parts), 0, 1)
    return np.prod(np.array(parts, dtype=float), axis=0)[()]


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
    parts = np.broadcast_arrays(*(np.asarray(x) for x in r))
    _interval("r", np.array(parts), 0, 1)
    return (1 - np.prod(1 - np.array(parts, dtype=float), axis=0))[()]


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


def _check_network(components, source, sink):
    """Validate a block diagram given as ``{name: (node, node)}``."""
    if not isinstance(components, Mapping) or not components:
        raise ValueError("components must map component names to (node, node) pairs")
    for name, ends in components.items():
        if not isinstance(name, str) or not isinstance(ends, tuple) or len(ends) != 2:
            raise ValueError("components must map component names (str) to (node, node) pairs")
    if source == sink:
        raise ValueError("source and sink must be different nodes")


def minimal_path_sets(components, source, sink):
    r"""Minimal path sets of a reliability block diagram drawn as a network.

    Each component is a link between two nodes, usable in either direction;
    the system works if the working components connect ``source`` to
    ``sink``. A minimal path set is a set of components whose working alone
    connects them and none of which can be dropped: the components along
    each simple path.

    Parameters
    ----------
    components : mapping of str to tuple
        ``{name: (node_a, node_b)}``; nodes are any hashable labels.
    source, sink : hashable
        Input and output nodes of the diagram.

    Returns
    -------
    list of tuple of str
        Minimal path sets, each sorted, shortest first.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    J. D. Esary and F. Proschan, "Coherent structures of non-identical
    components," *Technometrics* 5 (1963) 191-209.

    Examples
    --------
    >>> bridge = {"A": ("s", 1), "B": ("s", 2), "C": (1, 2), "D": (1, "t"), "E": (2, "t")}
    >>> minimal_path_sets(bridge, "s", "t")
    [('A', 'D'), ('B', 'E'), ('A', 'C', 'E'), ('B', 'C', 'D')]
    """
    _check_network(components, source, sink)
    links: dict = {}
    for name, (a, b) in components.items():
        links.setdefault(a, []).append((name, b))
        links.setdefault(b, []).append((name, a))
    found = set()

    def walk(node, visited, used):
        if node == sink:
            found.add(tuple(sorted(used)))
            return
        for name, other in links.get(node, []):
            if other not in visited:
                walk(other, visited | {other}, used + [name])

    walk(source, {source}, [])
    return sorted(found, key=lambda path: (len(path), path))


def minimal_cut_sets(components, source, sink):
    r"""Minimal cut sets of a reliability block diagram drawn as a network.

    A cut set is a set of components whose failure alone disconnects
    ``source`` from ``sink``; it is minimal if no component can be dropped.
    The minimal cut sets are the minimal sets that meet every minimal path
    set, found here by trying sets in order of size.

    Parameters
    ----------
    components : mapping of str to tuple
        ``{name: (node_a, node_b)}``.
    source, sink : hashable
        Input and output nodes of the diagram.

    Returns
    -------
    list of tuple of str
        Minimal cut sets, each sorted, smallest first. Empty if ``source``
        and ``sink`` are not connected at all.

    Notes
    -----
    Limited to 20 components that lie on some path. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    J. D. Esary and F. Proschan, "Coherent structures of non-identical
    components," *Technometrics* 5 (1963) 191-209.

    Examples
    --------
    >>> bridge = {"A": ("s", 1), "B": ("s", 2), "C": (1, 2), "D": (1, "t"), "E": (2, "t")}
    >>> minimal_cut_sets(bridge, "s", "t")
    [('A', 'B'), ('D', 'E'), ('A', 'C', 'E'), ('B', 'C', 'D')]
    """
    paths = [set(path) for path in minimal_path_sets(components, source, sink)]
    names = sorted(set().union(*paths))
    if len(names) > _MAX_COMPONENTS:
        raise ValueError(f"cut-set enumeration is limited to {_MAX_COMPONENTS} components")
    cuts: list[set] = []
    for size in range(1, len(names) + 1):
        for candidate in itertools.combinations(names, size):
            chosen = set(candidate)
            if all(chosen & path for path in paths) and not any(cut <= chosen for cut in cuts):
                cuts.append(chosen)
    return [tuple(sorted(cut)) for cut in cuts]


def _component_reliabilities(components, reliabilities):
    if not isinstance(reliabilities, Mapping):
        raise ValueError("reliabilities must map component names to reliabilities")
    missing = sorted(set(components) - set(reliabilities))
    if missing:
        raise ValueError(f"no reliability given for components {missing}")
    r = {name: float(reliabilities[name]) for name in components}
    _interval("reliabilities", list(r.values()), 0, 1)
    return r


def network_reliability(components, reliabilities, source, sink):
    r"""Exact reliability of a block diagram with independent components.

    The system works in every state of the components in which some minimal
    path set works entirely; the probabilities of those states are summed
    over all :math:`2^m` states of the :math:`m` components on paths.

    Parameters
    ----------
    components : mapping of str to tuple
        ``{name: (node_a, node_b)}``.
    reliabilities : mapping of str to float
        Reliability of each component [-], in [0, 1].
    source, sink : hashable
        Input and output nodes of the diagram.

    Returns
    -------
    float
        System reliability [-].

    Notes
    -----
    Limited to 20 components that lie on some path. Series, parallel and
    k-out-of-n systems are special cases. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/reliability`.

    Examples
    --------
    >>> bridge = {"A": ("s", 1), "B": ("s", 2), "C": (1, 2), "D": (1, "t"), "E": (2, "t")}
    >>> round(network_reliability(bridge, dict.fromkeys(bridge, 0.9), "s", "t"), 5)
    0.97848
    """
    r = _component_reliabilities(components, reliabilities)
    paths = minimal_path_sets(components, source, sink)
    names = sorted({name for path in paths for name in path})
    if len(names) > _MAX_COMPONENTS:
        raise ValueError(f"exact evaluation is limited to {_MAX_COMPONENTS} components")
    if not paths:
        return 0.0
    index = {name: i for i, name in enumerate(names)}
    codes = np.arange(2 ** len(names))
    states = (codes[:, None] >> np.arange(len(names))) & 1 == 1
    p = np.array([r[name] for name in names])
    weights = np.prod(np.where(states, p, 1 - p), axis=1)
    works = np.logical_or.reduce([np.all(states[:, [index[c] for c in path]], axis=1) for path in paths])
    return float(weights[works].sum())


def esary_proschan_bounds(components, reliabilities, source, sink):
    r"""Esary-Proschan lower and upper bounds on block-diagram reliability from minimal cuts and paths.

    :math:`\prod_{C}\left(1 - \prod_{i \in C} q_i\right) \le R \le
    1 - \prod_{P}\left(1 - \prod_{i \in P} p_i\right)`, over the minimal cut
    sets :math:`C` and path sets :math:`P`. The bounds treat the cuts (or
    paths) as if they failed independently; they are tight when components
    are reliable (lower) or unreliable (upper).

    Parameters
    ----------
    components : mapping of str to tuple
        ``{name: (node_a, node_b)}``.
    reliabilities : mapping of str to float
        Reliability of each component [-], in [0, 1].
    source, sink : hashable
        Input and output nodes of the diagram.

    Returns
    -------
    tuple of float
        ``(lower, upper)`` [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    J. D. Esary and F. Proschan, "Coherent structures of non-identical
    components," *Technometrics* 5 (1963) 191-209.

    Examples
    --------
    >>> bridge = {"A": ("s", 1), "B": ("s", 2), "C": (1, 2), "D": (1, "t"), "E": (2, "t")}
    >>> [round(b, 4) for b in esary_proschan_bounds(bridge, dict.fromkeys(bridge, 0.9), "s", "t")]
    [0.9781, 0.9973]
    """
    r = _component_reliabilities(components, reliabilities)
    paths = minimal_path_sets(components, source, sink)
    cuts = minimal_cut_sets(components, source, sink)
    lower = float(np.prod([1 - np.prod([1 - r[c] for c in cut]) for cut in cuts])) if paths else 0.0
    upper = 1 - float(np.prod([1 - np.prod([r[c] for c in path]) for path in paths]))
    return lower, upper


def standby_reliability(t, lam, n=2, switch_reliability=1.0):
    r"""Reliability of ``n`` identical units in cold standby, :math:`e^{-\lambda t}\sum_{k=0}^{n-1} p^k(\lambda t)^k/k!`.

    One unit runs; when it fails a switch, which works with probability
    :math:`p` at each changeover, brings in the next. Spares do not age
    while waiting. With a perfect switch the system life is the sum of
    :math:`n` exponential lives, a gamma (Erlang) distribution of shape
    :math:`n`.

    Parameters
    ----------
    t : float or array_like
        Mission time [s], nonnegative.
    lam : float
        Failure rate of a running unit [1/s], nonnegative.
    n : int, optional
        Number of units, the running one included, at least one.
    switch_reliability : float, optional
        Probability :math:`p` that a changeover succeeds [-], in [0, 1].

    Returns
    -------
    float or ndarray
        System reliability [-].

    Notes
    -----
    Cold standby with a perfect switch beats active parallel redundancy of
    the same units, since the spares do not wear while waiting. Inputs are
    raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    C. E. Ebeling, *An Introduction to Reliability and Maintainability
    Engineering* (McGraw-Hill, 1997), §6.5.

    Examples
    --------
    >>> round(float(standby_reliability(1000.0, 1e-3, 2)), 4)
    0.7358
    """
    _finite(t=t, lam=lam)
    _nonnegative(t=t, lam=lam)
    _integer("n", n, 1)
    _interval("switch_reliability", switch_reliability, 0, 1)
    x = lam * np.asarray(t, dtype=float)
    terms = sum((switch_reliability * x) ** k / math.factorial(k) for k in range(n))
    return np.exp(-x) * terms


def standby_mttf(lam, n=2, switch_reliability=1.0):
    r"""Mean life of ``n`` identical units in cold standby, :math:`\sum_{k=0}^{n-1} p^k/\lambda`.

    Parameters
    ----------
    lam : float
        Failure rate of a running unit [1/s], positive.
    n : int, optional
        Number of units, at least one.
    switch_reliability : float, optional
        Probability that a changeover succeeds [-], in [0, 1].

    Returns
    -------
    float
        Mean time to failure [s]; :math:`n/\lambda` with a perfect switch.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> standby_mttf(0.5, 3)
    6.0
    """
    _positive(lam=lam)
    _integer("n", n, 1)
    _interval("switch_reliability", switch_reliability, 0, 1)
    return sum(switch_reliability**k for k in range(n)) / lam
