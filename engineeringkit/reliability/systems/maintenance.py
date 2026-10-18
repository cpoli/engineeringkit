r"""Repairable systems: availability, Markov models and preventive replacement.

A repairable unit alternates between working and being repaired, so the
useful question is the fraction of time it is up. Markov models answer it
for constant failure and repair rates, from the transient after start-up to
the steady state :math:`\mathrm{MTBF}/(\mathrm{MTBF} + \mathrm{MTTR})`.
Barlow and Proschan showed when replacing a unit before it fails pays off:
only when its hazard rises.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg, special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.reliability.systems.lifetimes import weibull_reliability

__all__ = [
    "availability",
    "unavailability",
    "maintainability",
    "repair_rate",
    "markov_availability",
    "markov_state_probabilities",
    "age_replacement_cost_rate",
    "renewal_interval_cost_rate",
]


def availability(mtbf, mttr):
    r"""Steady-state availability :math:`A = \mathrm{MTBF}/(\mathrm{MTBF} + \mathrm{MTTR})`.

    Parameters
    ----------
    mtbf : float or array_like
        Mean time between failures [s], positive.
    mttr : float or array_like
        Mean time to repair [s], nonnegative.

    Returns
    -------
    float or ndarray
        Long-run fraction of time the unit is up [-].

    Notes
    -----
    This is a long-run average, not the probability of surviving a mission.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> availability(100.0, 25.0)
    0.8
    """
    _finite(mtbf=mtbf, mttr=mttr)
    _positive(mtbf=mtbf)
    _nonnegative(mttr=mttr)
    return mtbf / (mtbf + mttr)


def unavailability(mtbf, mttr):
    r"""Steady-state unavailability :math:`1 - A = \mathrm{MTTR}/(\mathrm{MTBF} + \mathrm{MTTR})`.

    Parameters
    ----------
    mtbf : float or array_like
        Mean time between failures [s], positive.
    mttr : float or array_like
        Mean time to repair [s], nonnegative.

    Returns
    -------
    float or ndarray
        Long-run fraction of time the unit is down [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(unavailability(100.0, 25.0), 2)
    0.2
    """
    return 1 - availability(mtbf, mttr)


def maintainability(t, mu):
    r"""Probability :math:`M(t) = 1 - e^{-\mu t}` that a repair is finished within ``t``.

    Parameters
    ----------
    t : float or array_like
        Time since the repair started [s], nonnegative.
    mu : float or array_like
        Constant repair rate [1/s], nonnegative.

    Returns
    -------
    float or ndarray
        Probability [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(float(maintainability(2.0, 0.5)), 4)
    0.6321
    """
    _finite(t=t, mu=mu)
    _nonnegative(t=t, mu=mu)
    return 1 - np.exp(-mu * np.asarray(t))


def repair_rate(mttr):
    r"""Constant repair rate :math:`\mu = 1/\mathrm{MTTR}`.

    Parameters
    ----------
    mttr : float or array_like
        Mean time to repair [s], positive.

    Returns
    -------
    float or ndarray
        Repair rate [1/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> repair_rate(2.0)
    0.5
    """
    _finite(mttr=mttr)
    _positive(mttr=mttr)
    return 1 / mttr


def markov_availability(t, lam, mu):
    r"""Point availability of a repairable unit that starts working at ``t = 0``.

    The two-state Markov model gives
    :math:`A(t) = \frac{\mu}{\lambda+\mu} + \frac{\lambda}{\lambda+\mu}e^{-(\lambda+\mu)t}`,
    which decays from one to the steady state with time constant
    :math:`1/(\lambda+\mu)`.

    Parameters
    ----------
    t : float or array_like
        Time [s], nonnegative.
    lam : float
        Constant failure rate [1/s], nonnegative.
    mu : float
        Constant repair rate [1/s], nonnegative; ``lam + mu`` must be positive.

    Returns
    -------
    float or ndarray
        Probability that the unit is up at time ``t`` [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    R. E. Barlow and F. Proschan, *Mathematical Theory of Reliability*
    (Wiley, 1965), ch. 5.

    Examples
    --------
    >>> round(float(markov_availability(1e6, 1e-3, 0.1)), 4)
    0.9901
    """
    _nonnegative(t=t, lam=lam, mu=mu)
    total = lam + mu
    _positive(rate_sum=total)
    return mu / total + lam / total * np.exp(-total * np.asarray(t))


def markov_state_probabilities(generator, initial, t):
    r"""State probabilities :math:`p(t) = p(0)\,e^{Qt}` of a continuous-time Markov chain.

    Parameters
    ----------
    generator : array_like, shape (n, n)
        Transition-rate matrix :math:`Q` [1/s]: nonnegative off-diagonal
        rates, rows summing to zero.
    initial : array_like, shape (n,)
        Initial state probabilities [-], nonnegative and summing to one.
    t : float or array_like
        Time(s) [s], nonnegative.

    Returns
    -------
    ndarray
        Shape ``(n,)`` for scalar ``t``, otherwise ``(len(t), n)`` [-].

    Notes
    -----
    Uses SciPy's matrix exponential. Availability is the summed probability
    of the up states. Inputs are raw SI numbers with no dimensional metadata;
    see :doc:`/guides/reliability`.

    Examples
    --------
    >>> Q = [[-1e-3, 1e-3], [0.1, -0.1]]
    >>> np.round(markov_state_probabilities(Q, [1, 0], 1e6), 4)
    array([0.9901, 0.0099])
    """
    Q = np.asarray(generator, dtype=float)
    if Q.ndim != 2 or Q.shape[0] != Q.shape[1] or not Q.size:
        raise ValueError("generator must be a nonempty square matrix")
    _finite(generator=Q)
    off = Q - np.diag(np.diag(Q))
    scale = float(np.abs(Q).max()) or 1.0
    if np.any(off < 0) or not np.allclose(Q.sum(axis=1), 0, atol=1e-12 * scale):
        raise ValueError("generator needs nonnegative off-diagonal rates and zero row sums")
    p0 = _vector("initial", initial)
    if p0.size != Q.shape[0]:
        raise ValueError("initial must have one probability per state")
    _nonnegative(initial=p0)
    if not np.isclose(p0.sum(), 1):
        raise ValueError("initial probabilities must sum to one")
    _nonnegative(t=t)
    times = np.atleast_1d(np.asarray(t, dtype=float))
    p = np.array([p0 @ linalg.expm(Q * ti) for ti in times])
    return p[0] if np.ndim(t) == 0 else p


def age_replacement_cost_rate(T, preventive_cost, failure_cost, eta, beta):
    r"""Long-run cost per unit time of replacing a Weibull unit at age ``T`` or at failure.

    Barlow and Proschan's renewal-reward ratio,
    :math:`C(T) = [c_p R(T) + c_f F(T)] / \int_0^T R(t)\,dt`, with the Weibull
    integral :math:`\eta\,\Gamma(1+1/\beta)\,P(1/\beta, (T/\eta)^\beta)`
    (:math:`P` the regularised lower incomplete gamma function).

    Parameters
    ----------
    T : float or array_like
        Planned replacement age [s], positive.
    preventive_cost : float
        Cost of a planned replacement [currency], nonnegative.
    failure_cost : float
        Cost of a replacement after failure [currency], nonnegative.
    eta : float
        Weibull scale [s], positive.
    beta : float
        Weibull shape [-], positive.

    Returns
    -------
    float or ndarray
        Cost rate [currency/s].

    Notes
    -----
    A finite optimum exists only if :math:`\beta > 1` and
    :math:`c_f > c_p`; otherwise the cost rate keeps falling toward
    :math:`c_f/\mathrm{MTTF}` (run to failure). Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    R. E. Barlow and F. Proschan, *Mathematical Theory of Reliability*
    (Wiley, 1965), ch. 4.

    Examples
    --------
    >>> round(float(age_replacement_cost_rate(1e9, 1.0, 10.0, 1000.0, 1.0)), 4)
    0.01
    """
    _positive(T=T, eta=eta, beta=beta)
    _nonnegative(preventive_cost=preventive_cost, failure_cost=failure_cost)
    R = weibull_reliability(T, eta, beta)
    uptime = eta * special.gamma(1 + 1 / beta) * special.gammainc(1 / beta, (np.asarray(T) / eta) ** beta)
    return (preventive_cost * R + failure_cost * (1 - R)) / uptime


def renewal_interval_cost_rate(preventive_cost, failure_cost, lam, T):
    r"""Constant-hazard maintenance cost-rate baseline :math:`c_p/T + c_f\lambda`.

    Parameters
    ----------
    preventive_cost : float
        Cost of each scheduled maintenance action [currency], nonnegative.
    failure_cost : float
        Cost of each failure [currency], nonnegative.
    lam : float
        Constant failure rate [1/s], nonnegative.
    T : float
        Maintenance interval [s], positive.

    Returns
    -------
    float
        Cost rate [currency/s].

    Notes
    -----
    With a constant hazard, scheduled replacement does not change the
    failure rate, so this cost only falls as ``T`` grows. Use
    :func:`age_replacement_cost_rate` for wear-out. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/reliability`.

    Examples
    --------
    >>> renewal_interval_cost_rate(10.0, 20.0, 0.1, 2.0)
    7.0
    """
    _finite(preventive_cost=preventive_cost, failure_cost=failure_cost, lam=lam, T=T)
    _nonnegative(preventive_cost=preventive_cost, failure_cost=failure_cost, lam=lam)
    _positive(T=T)
    return preventive_cost / T + failure_cost * lam
