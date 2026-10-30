r"""Probabilistic risk assessment: FMEA, fault trees, event trees and safety integrity levels.

MIL-P-1629 (1949) set out failure mode, effects and criticality analysis:
list every way each part can fail, its effect and how bad it is, and rank
the modes. Watson's fault tree works down from an undesired top event through AND and
OR gates to basic component failures. The Reactor Safety Study (WASH-1400)
worked forward from an initiating event through the success or failure of
each safety function. IEC 61508 grades a safety function by its average
probability of failure on demand.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.reliability.core.base import EventTreeResult, FMEAResult

__all__ = [
    "risk_priority_number",
    "fmea_ranking",
    "failure_mode_criticality",
    "fault_tree_probability",
    "event_tree_probability",
    "safety_integrity_pfd",
    "safety_integrity_level",
]

_MAX_ENUMERATED_EVENTS = 20


def _basic_events(tree, found):
    """Append every basic-event name in ``tree`` to ``found`` (with repeats) and check the gates."""
    if isinstance(tree, str):
        found.append(tree)
        return found
    if not isinstance(tree, tuple) or len(tree) < 2 or tree[0] not in ("and", "or"):
        raise ValueError("a fault tree is a basic-event name or a tuple ('and' | 'or', child, child, ...)")
    for child in tree[1:]:
        _basic_events(child, found)
    return found


def _gate_probability(tree, p):
    """Gate-by-gate probability, exact when no basic event appears twice."""
    if isinstance(tree, str):
        return p[tree]
    children = [_gate_probability(child, p) for child in tree[1:]]
    if tree[0] == "and":
        return float(np.prod(children))
    return 1 - float(np.prod([1 - c for c in children]))


def _occurs(tree, states, index):
    """Top-event indicator for every row of a boolean basic-event state table."""
    if isinstance(tree, str):
        return states[:, index[tree]]
    children = [_occurs(child, states, index) for child in tree[1:]]
    return np.logical_and.reduce(children) if tree[0] == "and" else np.logical_or.reduce(children)


def _ranking(name, value):
    """Check FMEA rankings: integers from 1 to 10, scalar or array."""
    a = np.asarray(value)
    _finite(**{name: a})
    if not np.issubdtype(a.dtype, np.integer) or np.any((a < 1) | (a > 10)):
        raise ValueError(f"{name} must be integer rankings from 1 to 10")
    return a


def risk_priority_number(severity, occurrence, detection):
    r"""FMEA risk priority number :math:`\mathrm{RPN} = S \times O \times D`.

    Parameters
    ----------
    severity : int or array_like of int
        Severity of the failure effect, 1 (none) to 10 (hazardous without
        warning) [-].
    occurrence : int or array_like of int
        Likelihood of the failure mode, 1 (remote) to 10 (almost certain)
        [-].
    detection : int or array_like of int
        Chance that current controls miss the mode before it reaches the
        user, 1 (certain detection) to 10 (no chance) [-].

    Returns
    -------
    int or ndarray of int
        Risk priority number, 1 to 1000 [-].

    Notes
    -----
    The rankings are ordinal, so the product orders modes but its ratios
    mean little: a mode with :math:`S = 10` can have a lower RPN than a
    trivial but frequent one, which is why high severities are reviewed
    whatever their RPN. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/reliability`.

    References
    ----------
    P. D. T. O'Connor and A. Kleyner, *Practical Reliability Engineering*
    (5th ed., Wiley, 2012), §7.3.

    Examples
    --------
    >>> risk_priority_number(8, 3, 6)
    144
    """
    s, o, d = (_ranking(n, v) for n, v in (("severity", severity), ("occurrence", occurrence), ("detection", detection)))
    rpn = s * o * d
    return int(rpn) if np.ndim(rpn) == 0 else rpn


def fmea_ranking(modes):
    r"""Rank failure modes by risk priority number, highest first.

    Parameters
    ----------
    modes : mapping of str to tuple of int
        ``{mode: (severity, occurrence, detection)}``, each ranking from 1
        to 10.

    Returns
    -------
    FMEAResult
        The modes and their rankings sorted by decreasing RPN (ties broken
        by severity, then name), with ``plot()`` for a Pareto chart.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    U.S. Department of Defense, *Procedures for Performing a Failure Mode,
    Effects and Criticality Analysis* (MIL-P-1629, 9 November 1949).

    Examples
    --------
    >>> result = fmea_ranking({"seal leak": (6, 5, 4), "bearing seizure": (9, 2, 7)})
    >>> result.modes, result.rpn.tolist()
    (('bearing seizure', 'seal leak'), [126, 120])
    """
    if not isinstance(modes, Mapping) or not modes:
        raise ValueError("modes must map failure-mode names to (severity, occurrence, detection)")
    rows = []
    for name, ranks in modes.items():
        if not isinstance(name, str) or not isinstance(ranks, tuple) or len(ranks) != 3:
            raise ValueError("modes must map failure-mode names (str) to (severity, occurrence, detection)")
        rows.append((risk_priority_number(*ranks), *ranks, name))
    rows.sort(key=lambda row: (-row[0], -row[1], row[4]))
    rpn, severity, occurrence, detection, names = zip(*rows)
    return FMEAResult(tuple(names), np.array(severity), np.array(occurrence), np.array(detection), np.array(rpn))


def failure_mode_criticality(mode_ratio, loss_probability, lam, time):
    r"""MIL-STD-1629A failure-mode criticality number :math:`C_m = \beta\,\alpha\,\lambda_p\,t`.

    Parameters
    ----------
    mode_ratio : float or array_like
        Fraction :math:`\alpha` of the part's failures that occur in this
        mode [-], in [0, 1].
    loss_probability : float or array_like
        Conditional probability :math:`\beta` that the mode causes the
        effect being ranked (loss of mission) [-], in [0, 1].
    lam : float or array_like
        Part failure rate :math:`\lambda_p` [1/s], nonnegative.
    time : float or array_like
        Operating time per mission :math:`t` [s], nonnegative.

    Returns
    -------
    float or ndarray
        Expected number of losses from this mode per mission [-]; summing
        over a part's modes gives the item criticality.

    Notes
    -----
    This is the quantitative criticality of the military standard; the
    RPN of :func:`risk_priority_number` is the qualitative ranking used
    where failure rates are unknown. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    U.S. Department of Defense, *Procedures for Performing a Failure Mode,
    Effects and Criticality Analysis* (MIL-STD-1629A, 1980), task 102.

    Examples
    --------
    >>> round(float(failure_mode_criticality(0.4, 0.5, 2e-6, 1e4)), 6)
    0.004
    """
    _interval("mode_ratio", mode_ratio, 0, 1)
    _interval("loss_probability", loss_probability, 0, 1)
    _nonnegative(lam=lam, time=time)
    return np.asarray(loss_probability) * mode_ratio * lam * time


def fault_tree_probability(tree, probabilities):
    r"""Probability of a fault tree's top event from independent basic-event probabilities.

    A tree is a basic-event name (``str``) or a gate tuple
    ``("and", child, ...)`` / ``("or", child, ...)``. If no basic event
    appears twice, gates combine directly: :math:`\prod p_i` for AND and
    :math:`1 - \prod(1 - p_i)` for OR. A repeated event makes the gate inputs
    dependent, so the top event is then evaluated exactly over every state
    of the basic events.

    Parameters
    ----------
    tree : str or tuple
        Fault-tree structure.
    probabilities : mapping of str to float
        Probability of each basic event [-], in [0, 1].

    Returns
    -------
    float
        Top-event probability [-].

    Notes
    -----
    Exact enumeration is limited to 20 distinct basic events. Inputs are raw
    SI numbers with no dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    H. A. Watson, *Launch Control Safety Study* (Bell Telephone
    Laboratories, 1961); D. F. Haasl, "Advanced concepts in fault tree
    analysis," *System Safety Symposium* (Seattle, 1965).

    Examples
    --------
    >>> round(fault_tree_probability(("or", "A", ("and", "B", "C")), {"A": 0.1, "B": 0.2, "C": 0.5}), 3)
    0.19
    """
    events = _basic_events(tree, [])
    if not isinstance(probabilities, Mapping):
        raise ValueError("probabilities must map basic-event names to probabilities")
    missing = sorted(set(events) - set(probabilities))
    if missing:
        raise ValueError(f"no probability given for basic events {missing}")
    names = sorted(set(events))
    p = {name: float(probabilities[name]) for name in names}
    _interval("probabilities", list(p.values()), 0, 1)
    if len(names) == len(events):
        return _gate_probability(tree, p)
    if len(names) > _MAX_ENUMERATED_EVENTS:
        raise ValueError(f"exact evaluation of repeated events is limited to {_MAX_ENUMERATED_EVENTS} basic events")
    codes = np.arange(2 ** len(names))
    states = (codes[:, None] >> np.arange(len(names))) & 1 == 1
    q = np.array([p[name] for name in names])
    weights = np.prod(np.where(states, q, 1 - q), axis=1)
    return float(weights[_occurs(tree, states, {name: i for i, name in enumerate(names)})].sum())


def event_tree_probability(initiating_frequency, failure_probabilities):
    r"""Frequencies of every sequence of an event tree with independent branch points.

    After an initiating event of frequency :math:`f_0`, each of :math:`m`
    safety functions succeeds or fails; a sequence's frequency is
    :math:`f_0 \prod_j q_j^{\,o_j}(1-q_j)^{1-o_j}`.

    Parameters
    ----------
    initiating_frequency : float
        Initiating-event frequency [1/s] (or a probability), nonnegative.
    failure_probabilities : array_like
        Failure-on-demand probability of each safety function [-], in
        [0, 1], in the order the tree asks them (at most 20).

    Returns
    -------
    EventTreeResult
        The :math:`2^m` outcome rows (first function most significant, all
        successes first) and their frequencies.

    Notes
    -----
    Real event trees prune branches that cannot matter and condition later
    functions on earlier ones; this helper keeps every branch and assumes
    independence. Inputs are raw SI numbers with no dimensional metadata;
    see :doc:`/guides/reliability`.

    References
    ----------
    U.S. Nuclear Regulatory Commission, *Reactor Safety Study* (WASH-1400,
    NUREG-75/014, 1975).

    Examples
    --------
    >>> result = event_tree_probability(1e-2, [0.1, 0.5])
    >>> np.round(result.frequencies, 4)
    array([0.0045, 0.0045, 0.0005, 0.0005])
    """
    _nonnegative(initiating_frequency=initiating_frequency)
    q = _vector("failure_probabilities", failure_probabilities)
    _interval("failure_probabilities", q, 0, 1)
    if q.size > _MAX_ENUMERATED_EVENTS:
        raise ValueError(f"at most {_MAX_ENUMERATED_EVENTS} safety functions")
    codes = np.arange(2**q.size)
    outcomes = (codes[:, None] >> np.arange(q.size)[::-1]) & 1 == 1
    frequencies = initiating_frequency * np.prod(np.where(outcomes, q, 1 - q), axis=1)
    return EventTreeResult(outcomes, frequencies)


def safety_integrity_pfd(lam_du, test_interval):
    r"""Average probability of failure on demand of a single channel, :math:`\lambda_{DU} T/2`.

    Parameters
    ----------
    lam_du : float or array_like
        Dangerous undetected failure rate [1/s], nonnegative.
    test_interval : float or array_like
        Proof-test interval [s], positive.

    Returns
    -------
    float or ndarray
        PFDavg [-].

    Raises
    ------
    ValueError
        If any :math:`\lambda_{DU} T > 0.1`, where the rare-failure
        approximation no longer holds.

    Notes
    -----
    Assumes perfect proof tests and negligible repair time. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    IEC 61508-6, *Functional safety of electrical/electronic/programmable
    electronic safety-related systems* (1998; 2nd ed. 2010), Annex B.

    Examples
    --------
    >>> safety_integrity_pfd(0.001, 10.0)
    0.005
    """
    _finite(lam_du=lam_du, test_interval=test_interval)
    _nonnegative(lam_du=lam_du)
    _positive(test_interval=test_interval)
    product = np.asarray(lam_du) * np.asarray(test_interval)
    if np.any(product > 0.1):
        raise ValueError("rare-failure approximation requires lambda*T <= 0.1")
    return product / 2 if product.ndim else float(product) / 2


def safety_integrity_level(pfd_avg):
    r"""IEC 61508 safety integrity level for a low-demand PFDavg (0 if it meets none).

    ============  =====================
    SIL           PFDavg
    ============  =====================
    4             :math:`< 10^{-4}`
    3             :math:`[10^{-4}, 10^{-3})`
    2             :math:`[10^{-3}, 10^{-2})`
    1             :math:`[10^{-2}, 10^{-1})`
    0             :math:`\ge 10^{-1}`
    ============  =====================

    Parameters
    ----------
    pfd_avg : float or array_like
        Average probability of failure on demand [-], in [0, 1].

    Returns
    -------
    int or ndarray of int
        Safety integrity level.

    Notes
    -----
    The standard's SIL 4 band starts at :math:`10^{-5}`; a lower PFDavg is
    still reported as SIL 4. A level also needs architectural constraints
    and systematic capability, which a number alone cannot show. Inputs are
    raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> safety_integrity_level(8.76e-3)
    2
    """
    _interval("pfd_avg", pfd_avg, 0, 1)
    level = 4 - np.digitize(pfd_avg, [1e-4, 1e-3, 1e-2, 1e-1])
    return int(level) if np.ndim(level) == 0 else level
