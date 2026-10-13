"""Reliability engineering: lifetime distributions, systems, availability and maintenance."""

import math

import numpy as np
from scipy import special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive


def exponential_reliability(t, lam):
    """Exponential reliability.

    Parameters
    ----------
    t : s
        Input in the units stated above.
    lam : 1/s
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(t=t, lam=lam)
    _nonnegative(t=t, lam=lam)
    return np.exp(-lam * np.asarray(t))


def exponential_failure_probability(t, lam):
    """Exponential failure probability.

    Parameters
    ----------
    t : s
        Input in the units stated above.
    lam : 1/s
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _nonnegative(t=t, lam=lam)
    return -np.expm1(-lam * np.asarray(t))


def exponential_mttf(lam):
    """Exponential mttf.

    Parameters
    ----------
    lam : 1/s
        Input in the units stated above.

    Returns
    -------
    result
        s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(lam=lam)
    _positive(lam=lam)
    return 1 / lam


def exponential_hazard(lam):
    """Exponential hazard.

    Parameters
    ----------
    lam : 1/s
        Input in the units stated above.

    Returns
    -------
    result
        1/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(lam=lam)
    _nonnegative(lam=lam)
    return lam


def weibull_reliability(t, eta, beta):
    """Weibull reliability.

    Parameters
    ----------
    t : s
        Input in the units stated above.
    eta : s
        Input in the units stated above.
    beta : 1
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(t=t, eta=eta, beta=beta)
    _nonnegative(t=t)
    _positive(eta=eta, beta=beta)
    return np.exp(-((np.asarray(t) / eta) ** beta))


def weibull_cdf(t, eta, beta):
    """Weibull cdf.

    Parameters
    ----------
    t : s
        Input in the units stated above.
    eta : s
        Input in the units stated above.
    beta : 1
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _nonnegative(t=t)
    _positive(eta=eta, beta=beta)
    return -np.expm1(-((np.asarray(t) / eta) ** beta))


def weibull_hazard(t, eta, beta):
    """Weibull hazard.

    Parameters
    ----------
    t : s
        Input in the units stated above.
    eta : s
        Input in the units stated above.
    beta : 1
        Input in the units stated above.

    Returns
    -------
    result
        1/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(t=t, eta=eta, beta=beta)
    _nonnegative(t=t)
    _positive(eta=eta, beta=beta)
    with np.errstate(divide="ignore"):
        return beta / eta * (np.asarray(t) / eta) ** (beta - 1)


def weibull_mttf(eta, beta):
    """Weibull mttf.

    Parameters
    ----------
    eta : s
        Input in the units stated above.
    beta : 1
        Input in the units stated above.

    Returns
    -------
    result
        s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(eta=eta, beta=beta)
    _positive(eta=eta, beta=beta)
    return eta * math.gamma(1 + 1 / beta)


def series_reliability(*r):
    """Series reliability.

    Parameters
    ----------
    r : 1
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _interval("r", r, 0, 1)
    return float(np.prod(r))


def parallel_reliability(*r):
    """Parallel reliability.

    Parameters
    ----------
    r : 1
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _interval("r", r, 0, 1)
    return 1 - float(np.prod([1 - x for x in r]))


def k_of_n_reliability(n, k, r):
    """K of n reliability.

    Parameters
    ----------
    n : 1
        Input in the units stated above.
    k : 1
        Input in the units stated above.
    r : 1
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(n=n, k=k, r=r)
    _integer("n", n, 1)
    _integer("k", k, 0)
    _interval("r", r, 0, 1)
    if k > n:
        raise ValueError("k cannot exceed n")
    return sum(math.comb(n, i) * r**i * (1 - r) ** (n - i) for i in range(k, n + 1))


def availability(mtbf, mttr):
    """Availability.

    Parameters
    ----------
    mtbf : s
        Input in the units stated above.
    mttr : s
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(mtbf=mtbf, mttr=mttr)
    _positive(mtbf=mtbf)
    _nonnegative(mttr=mttr)
    return mtbf / (mtbf + mttr)


def unavailability(mtbf, mttr):
    """Unavailability.

    Parameters
    ----------
    mtbf : s
        Input in the units stated above.
    mttr : s
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return 1 - availability(mtbf, mttr)


def failure_rate_from_mtbf(mtbf):
    """Failure rate from mtbf.

    Parameters
    ----------
    mtbf : s
        Input in the units stated above.

    Returns
    -------
    result
        1/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(mtbf=mtbf)
    _positive(mtbf=mtbf)
    return 1 / mtbf


def maintainability(t, mu):
    """Maintainability.

    Parameters
    ----------
    t : s
        Input in the units stated above.
    mu : 1/s
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(t=t, mu=mu)
    _nonnegative(t=t, mu=mu)
    return 1 - np.exp(-mu * np.asarray(t))


def repair_rate(mttr):
    """Repair rate.

    Parameters
    ----------
    mttr : s
        Input in the units stated above.

    Returns
    -------
    result
        1/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(mttr=mttr)
    _positive(mttr=mttr)
    return 1 / mttr


def expected_failures(lam, total_time):
    """Expected failures.

    Parameters
    ----------
    lam : 1/s
        Input in the units stated above.
    total_time : s
        Input in the units stated above.

    Returns
    -------
    result
        expected count.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(lam=lam, total_time=total_time)
    _nonnegative(lam=lam, total_time=total_time)
    return lam * total_time


def poisson_failures_probability(k, lam, total_time):
    """Evaluate a Poisson count probability stably in log space.

    Parameters
    ----------
    k : 1
        Input in the units stated above.
    lam : 1/s
        Input in the units stated above.
    total_time : s
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(k=k, lam=lam, total_time=total_time)
    _integer("k", k, 0)
    _nonnegative(lam=lam, total_time=total_time)
    mean = lam * total_time
    return float(np.exp(special.xlogy(k, mean) - mean - special.gammaln(k + 1)))


def birnbaum_importance_series(reliabilities, i):
    """Birnbaum importance series.

    Parameters
    ----------
    reliabilities : 1
        Input in the units stated above.
    i : 1
        Input in the units stated above.

    Returns
    -------
    result
        dimensionless sensitivity.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(reliabilities=reliabilities, i=i)
    _interval("reliabilities", reliabilities, 0, 1)
    _integer("i", i, 0)
    if i >= len(reliabilities):
        raise ValueError("component index out of range")
    return float(np.prod(np.delete(np.asarray(reliabilities), i)))


def reliability_growth_duane(t, K, alpha):
    """Reliability growth duane.

    Parameters
    ----------
    t : s
        Input in the units stated above.
    K : s^(1-alpha)
        Input in the units stated above.
    alpha : 1
        Input in the units stated above.

    Returns
    -------
    result
        s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(t=t, K=K, alpha=alpha)
    _positive(t=t, K=K)
    _nonnegative(alpha=alpha)
    return K * np.asarray(t) ** alpha


def bathtub_hazard(t, early_scale=1, early_rate=1, random_rate=0.01, wear_scale=1e-6, wear_power=2):
    """Bathtub hazard.

    Parameters
    ----------
    t : s
        Input in the units stated above.
    early_scale : 1/s
        Input in the units stated above.
    early_rate : 1/s
        Input in the units stated above.
    random_rate : 1/s
        Input in the units stated above.
    wear_scale : 1/s^(wear_power+1)
        Input in the units stated above.
    wear_power : 1
        Input in the units stated above.

    Returns
    -------
    result
        1/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(t=t, early_scale=early_scale, early_rate=early_rate, random_rate=random_rate, wear_scale=wear_scale, wear_power=wear_power)
    _nonnegative(t=t, early_scale=early_scale, early_rate=early_rate, random_rate=random_rate, wear_scale=wear_scale, wear_power=wear_power)
    t = np.asarray(t)
    return early_scale * np.exp(-early_rate * t) + random_rate + wear_scale * t**wear_power


def safety_integrity_pfd(lam_du, test_interval):
    """Estimate average dangerous undetected failure probability for rare failures.

    Parameters
    ----------
    lam_du : 1/s
        Input in the units stated above.
    test_interval : s
        Input in the units stated above.

    Returns
    -------
    result
        probability [0,1].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(lam_du=lam_du, test_interval=test_interval)
    _nonnegative(lam_du=lam_du)
    _positive(test_interval=test_interval)
    if lam_du * test_interval > 0.1:
        raise ValueError("rare-failure approximation requires lambda*T <= 0.1")
    return lam_du * test_interval / 2


def renewal_interval_cost_rate(preventive_cost, failure_cost, lam, T):
    """Compute a simple scheduled-maintenance cost-rate baseline, not a renewal solver.

    Parameters
    ----------
    preventive_cost : currency units
        Input in the units stated above.
    failure_cost : currency units
        Input in the units stated above.
    lam : 1/s
        Input in the units stated above.
    T : s
        Input in the units stated above.

    Returns
    -------
    result
        currency/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/reliability`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(preventive_cost=preventive_cost, failure_cost=failure_cost, lam=lam, T=T)
    _nonnegative(preventive_cost=preventive_cost, failure_cost=failure_cost, lam=lam)
    _positive(T=T)
    return preventive_cost / T + failure_cost * lam
