r"""Lifetime distributions, failure counts and life-test estimates.

Halley's life table turned survivor counts into age-specific death rates;
Poisson counted rare events; Weibull gave a two-parameter law flexible
enough for fatigue, bearings and electronics; Epstein and Sobel put
confidence limits on exponential life tests; and Kaplan and Meier estimated
a survival curve from tests in which some units never failed. Times are in
seconds and rates in 1/s, but any consistent unit (hours) works.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import special, stats

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.reliability.core.base import KaplanMeierResult

__all__ = [
    "exponential_reliability",
    "exponential_failure_probability",
    "exponential_mttf",
    "exponential_hazard",
    "failure_rate_from_mtbf",
    "weibull_reliability",
    "weibull_cdf",
    "weibull_hazard",
    "weibull_mttf",
    "bathtub_hazard",
    "life_table_hazard",
    "expected_failures",
    "poisson_failures_probability",
    "mtbf_confidence_bounds",
    "kaplan_meier",
]


def exponential_reliability(t, lam):
    r"""Survival probability :math:`R(t) = e^{-\lambda t}` under a constant hazard.

    Parameters
    ----------
    t : float or array_like
        Time [s], nonnegative.
    lam : float or array_like
        Constant failure rate :math:`\lambda` [1/s], nonnegative.

    Returns
    -------
    float or ndarray
        Reliability [-], in [0, 1].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(float(exponential_reliability(1000.0, 1e-3)), 4)
    0.3679
    """
    _finite(t=t, lam=lam)
    _nonnegative(t=t, lam=lam)
    return np.exp(-lam * np.asarray(t))


def exponential_failure_probability(t, lam):
    r"""Failure probability :math:`F(t) = 1 - e^{-\lambda t}`, accurate for tiny :math:`\lambda t`.

    Parameters
    ----------
    t : float or array_like
        Time [s], nonnegative.
    lam : float or array_like
        Constant failure rate [1/s], nonnegative.

    Returns
    -------
    float or ndarray
        Probability of failure by time ``t`` [-].

    Notes
    -----
    ``expm1`` keeps the result accurate when :math:`\lambda t \ll 1`.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(float(exponential_failure_probability(1.0, 1e-12)) / 1e-12, 9)
    1.0
    """
    _nonnegative(t=t, lam=lam)
    return -np.expm1(-lam * np.asarray(t))


def exponential_mttf(lam):
    r"""Mean time to failure :math:`1/\lambda` of an exponential life.

    Parameters
    ----------
    lam : float
        Constant failure rate [1/s], positive.

    Returns
    -------
    float
        Mean time to failure [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> exponential_mttf(0.5)
    2.0
    """
    _finite(lam=lam)
    _positive(lam=lam)
    return 1 / lam


def exponential_hazard(lam):
    r"""Hazard rate of an exponential life, the constant :math:`\lambda` itself.

    Parameters
    ----------
    lam : float or array_like
        Constant failure rate [1/s], nonnegative.

    Returns
    -------
    float or ndarray
        Hazard rate [1/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> exponential_hazard(2e-6)
    2e-06
    """
    _finite(lam=lam)
    _nonnegative(lam=lam)
    return lam


def failure_rate_from_mtbf(mtbf):
    r"""Constant failure rate :math:`\lambda = 1/\mathrm{MTBF}`.

    Parameters
    ----------
    mtbf : float or array_like
        Mean time between failures [s], positive.

    Returns
    -------
    float or ndarray
        Failure rate [1/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> failure_rate_from_mtbf(10.0)
    0.1
    """
    _finite(mtbf=mtbf)
    _positive(mtbf=mtbf)
    return 1 / mtbf


def weibull_reliability(t, eta, beta):
    r"""Weibull survival probability :math:`R(t) = e^{-(t/\eta)^\beta}`.

    Parameters
    ----------
    t : float or array_like
        Time [s], nonnegative.
    eta : float or array_like
        Scale (characteristic life, 63.2% failed) [s], positive.
    beta : float or array_like
        Shape [-], positive: below one the hazard falls, at one it is
        constant (exponential), above one it rises.

    Returns
    -------
    float or ndarray
        Reliability [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    W. Weibull, "A statistical distribution function of wide
    applicability," *J. Appl. Mech.* 18 (1951) 293-297.

    Examples
    --------
    >>> round(float(weibull_reliability(3.0, 3.0, 2.0)), 4)
    0.3679
    """
    _finite(t=t, eta=eta, beta=beta)
    _nonnegative(t=t)
    _positive(eta=eta, beta=beta)
    return np.exp(-((np.asarray(t) / eta) ** beta))


def weibull_cdf(t, eta, beta):
    r"""Weibull failure probability :math:`F(t) = 1 - e^{-(t/\eta)^\beta}`.

    Parameters
    ----------
    t : float or array_like
        Time [s], nonnegative.
    eta : float or array_like
        Scale [s], positive.
    beta : float or array_like
        Shape [-], positive.

    Returns
    -------
    float or ndarray
        Probability of failure by time ``t`` [-].

    Notes
    -----
    ``expm1`` keeps early-life probabilities accurate. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(float(weibull_cdf(1.0, 1.0, 1.0)), 4)
    0.6321
    """
    _nonnegative(t=t)
    _positive(eta=eta, beta=beta)
    return -np.expm1(-((np.asarray(t) / eta) ** beta))


def weibull_hazard(t, eta, beta):
    r"""Weibull hazard rate :math:`h(t) = (\beta/\eta)(t/\eta)^{\beta-1}`.

    Parameters
    ----------
    t : float or array_like
        Time [s], nonnegative.
    eta : float or array_like
        Scale [s], positive.
    beta : float or array_like
        Shape [-], positive.

    Returns
    -------
    float or ndarray
        Hazard rate [1/s]; infinite at ``t = 0`` when ``beta < 1``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> float(weibull_hazard(2.0, 2.0, 1.0))
    0.5
    """
    _finite(t=t, eta=eta, beta=beta)
    _nonnegative(t=t)
    _positive(eta=eta, beta=beta)
    with np.errstate(divide="ignore"):
        return beta / eta * (np.asarray(t) / eta) ** (beta - 1)


def weibull_mttf(eta, beta):
    r"""Weibull mean life :math:`\eta\,\Gamma(1 + 1/\beta)`.

    Parameters
    ----------
    eta : float or array_like
        Scale [s], positive.
    beta : float or array_like
        Shape [-], positive.

    Returns
    -------
    float or ndarray
        Mean time to failure [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(float(weibull_mttf(3.0, 2.0)), 5)
    2.65868
    """
    _finite(eta=eta, beta=beta)
    _positive(eta=eta, beta=beta)
    if np.ndim(eta) == 0 and np.ndim(beta) == 0:
        return eta * math.gamma(1 + 1 / beta)
    return np.asarray(eta) * special.gamma(1 + 1 / np.asarray(beta))


def bathtub_hazard(t, early_scale=1, early_rate=1, random_rate=0.01, wear_scale=1e-6, wear_power=2):
    r"""Bathtub-shaped hazard: decaying infant mortality, a constant floor and power-law wear-out.

    :math:`h(t) = a e^{-b t} + c + d\,t^{p}`.

    Parameters
    ----------
    t : float or array_like
        Age [s], nonnegative.
    early_scale : float, optional
        Infant-mortality hazard at age zero, :math:`a` [1/s].
    early_rate : float, optional
        Decay rate of infant mortality, :math:`b` [1/s].
    random_rate : float, optional
        Constant random-failure hazard, :math:`c` [1/s].
    wear_scale : float, optional
        Wear-out coefficient, :math:`d` [1/s^(p+1)].
    wear_power : float, optional
        Wear-out exponent, :math:`p` [-].

    Returns
    -------
    float or ndarray
        Hazard rate [1/s].

    Notes
    -----
    All parameters are nonnegative. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/reliability`.

    Examples
    --------
    >>> float(bathtub_hazard(0.0))
    1.01
    """
    params = dict(early_scale=early_scale, early_rate=early_rate, random_rate=random_rate, wear_scale=wear_scale, wear_power=wear_power)
    _finite(t=t, **params)
    _nonnegative(t=t, **params)
    t = np.asarray(t)
    return early_scale * np.exp(-early_rate * t) + random_rate + wear_scale * t**wear_power


def life_table_hazard(ages, survivors):
    r"""Piecewise-constant hazard rates from a life table of survivor counts.

    Over each interval :math:`[x_i, x_{i+1})` the hazard
    :math:`h_i = -\ln(l_{i+1}/l_i)/(x_{i+1}-x_i)` reproduces the recorded
    fall from :math:`l_i` to :math:`l_{i+1}` survivors exactly.

    Parameters
    ----------
    ages : array_like
        Ages at which survivors are counted [s], strictly increasing, at
        least two.
    survivors : array_like
        Number alive at each age [-], positive and nonincreasing.

    Returns
    -------
    ndarray
        Hazard rate on each of the ``len(ages) - 1`` intervals [1/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    E. Halley, "An estimate of the degrees of the mortality of mankind,"
    *Phil. Trans. R. Soc.* 17 (1693) 596-610.

    Examples
    --------
    >>> np.round(life_table_hazard([0, 1, 2], [1000, 500, 250]), 4)
    array([0.6931, 0.6931])
    """
    x = _vector("ages", ages, min_length=2)
    n = _vector("survivors", survivors, min_length=2)
    if x.shape != n.shape:
        raise ValueError("ages and survivors must have the same length")
    _positive(survivors=n)
    if np.any(np.diff(x) <= 0):
        raise ValueError("ages must be strictly increasing")
    if np.any(np.diff(n) > 0):
        raise ValueError("survivors must be nonincreasing")
    return -np.diff(np.log(n)) / np.diff(x)


def expected_failures(lam, total_time):
    r"""Expected number of failures :math:`\lambda T` in a total operating time.

    Parameters
    ----------
    lam : float or array_like
        Constant failure rate [1/s], nonnegative.
    total_time : float or array_like
        Accumulated operating time [s], nonnegative.

    Returns
    -------
    float or ndarray
        Expected count [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> expected_failures(0.5, 2.0)
    1.0
    """
    _finite(lam=lam, total_time=total_time)
    _nonnegative(lam=lam, total_time=total_time)
    return lam * total_time


def poisson_failures_probability(k, lam, total_time):
    r"""Poisson probability of exactly ``k`` failures, :math:`(\lambda T)^k e^{-\lambda T}/k!`.

    Parameters
    ----------
    k : int
        Number of failures, nonnegative.
    lam : float
        Constant failure rate [1/s], nonnegative.
    total_time : float
        Accumulated operating time [s], nonnegative.

    Returns
    -------
    float
        Probability [-].

    Notes
    -----
    Evaluated in log space with ``gammaln`` so that large counts do not
    overflow. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    S.-D. Poisson, *Recherches sur la probabilité des jugements en matière
    criminelle et en matière civile* (Paris, Bachelier, 1837).

    Examples
    --------
    >>> round(poisson_failures_probability(0, 1.0, 1.0), 4)
    0.3679
    """
    _finite(k=k, lam=lam, total_time=total_time)
    _integer("k", k, 0)
    _nonnegative(lam=lam, total_time=total_time)
    mean = lam * total_time
    return float(np.exp(special.xlogy(k, mean) - mean - special.gammaln(k + 1)))


def mtbf_confidence_bounds(total_time, failures, confidence=0.9, time_terminated=False):
    r"""Two-sided chi-square confidence bounds on an exponential MTBF from a life test.

    With :math:`r` failures in total unit-time :math:`T`, the estimate is
    :math:`\hat\theta = T/r` and :math:`2T/\theta` is chi-square with
    :math:`2r` degrees of freedom for a failure-terminated test (Epstein and
    Sobel). A time-terminated test uses :math:`2r + 2` degrees of freedom for
    the lower bound.

    Parameters
    ----------
    total_time : float
        Total time on test summed over all units, censored ones included [s],
        positive.
    failures : int
        Number of failures :math:`r`; at least one for a failure-terminated
        test, at least zero for a time-terminated one.
    confidence : float, optional
        Two-sided confidence level [-], in (0, 1).
    time_terminated : bool, optional
        ``True`` if the test stopped at a fixed time rather than at the
        :math:`r`-th failure.

    Returns
    -------
    tuple of float
        ``(lower, upper)`` MTBF bounds [s]; ``upper`` is infinite for a
        time-terminated test with no failures.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    B. Epstein and M. Sobel, "Life testing," *J. Amer. Statist. Assoc.* 48
    (1953) 486-502.

    Examples
    --------
    >>> [round(b) for b in mtbf_confidence_bounds(10000.0, 5)]
    [1092, 5076]
    """
    _positive(total_time=total_time)
    _integer("failures", failures, 0 if time_terminated else 1)
    _interval("confidence", confidence, 0, 1, closed=False)
    tail = (1 - confidence) / 2
    lower = 2 * total_time / stats.chi2.ppf(1 - tail, 2 * failures + (2 if time_terminated else 0))
    upper = 2 * total_time / stats.chi2.ppf(tail, 2 * failures) if failures else math.inf
    return float(lower), float(upper)


def kaplan_meier(times, failed=None):
    r"""Kaplan-Meier product-limit estimate of survival from censored life data.

    At each distinct failure time :math:`t_j`, with :math:`n_j` units at risk
    and :math:`d_j` failures, :math:`\hat S(t) = \prod_{t_j \le t}(1 - d_j/n_j)`.
    Censored units leave the risk set without moving the curve; a unit
    censored at a failure time counts as at risk for that failure.

    Parameters
    ----------
    times : array_like
        Failure or censoring time of each unit [s], nonnegative.
    failed : array_like of {0, 1} or bool, optional
        ``1`` where the unit failed, ``0`` where it was censored (removed or
        still running). Default: every unit failed.

    Returns
    -------
    KaplanMeierResult
        Failure times, survival estimate, risk-set sizes, failure counts and
        Greenwood standard errors.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    E. L. Kaplan and P. Meier, "Nonparametric estimation from incomplete
    observations," *J. Amer. Statist. Assoc.* 53 (1958) 457-481.

    Examples
    --------
    >>> result = kaplan_meier([1, 2, 3, 4], [1, 0, 1, 1])
    >>> np.round(result.survival, 4)
    array([0.75 , 0.375, 0.   ])
    """
    t = _vector("times", times)
    _nonnegative(times=t)
    if failed is None:
        d = np.ones(t.size, dtype=bool)
    else:
        flags = np.asarray(failed)
        if flags.shape != t.shape or not np.all((flags == 0) | (flags == 1)):
            raise ValueError("failed must be a 0/1 vector matching times")
        d = flags.astype(bool)
    u = np.unique(t[d])
    at_risk = np.array([np.count_nonzero(t >= x) for x in u], dtype=int)
    events = np.array([np.count_nonzero(d & (t == x)) for x in u], dtype=int)
    survival = np.cumprod(1 - events / at_risk)
    survivors = at_risk - events
    terms = np.divide(events, at_risk * survivors, out=np.zeros(u.size), where=survivors > 0)
    std_error = survival * np.sqrt(np.cumsum(terms))
    return KaplanMeierResult(u, survival, at_risk, events, std_error, int(t.size), int(np.count_nonzero(~d)))
