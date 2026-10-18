r"""Standard normal conversions between a reliability index and a failure probability.

Structural reliability states its answer as an index :math:`\beta`, the
number of standard deviations between the mean state and failure, with
:math:`p_f = \Phi(-\beta)`. SciPy's ``ndtr`` and ``ndtri`` evaluate the
normal distribution accurately far into the tail.
"""

from __future__ import annotations

import numpy as np
from scipy import special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval

__all__ = ["failure_probability_from_index", "reliability_index_from_probability"]


def failure_probability_from_index(beta):
    r"""Failure probability :math:`p_f = \Phi(-\beta)` for a reliability index.

    Parameters
    ----------
    beta : float or array_like
        Reliability index [-]; negative values mean failure is more likely
        than not.

    Returns
    -------
    float or ndarray
        Failure probability [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(float(failure_probability_from_index(3.0)), 6)
    0.00135
    """
    _finite(beta=beta)
    return special.ndtr(-np.asarray(beta, dtype=float))


def reliability_index_from_probability(pf):
    r"""Reliability index :math:`\beta = -\Phi^{-1}(p_f)` for a failure probability.

    Parameters
    ----------
    pf : float or array_like
        Failure probability [-], in the open interval (0, 1).

    Returns
    -------
    float or ndarray
        Reliability index [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    Examples
    --------
    >>> round(float(reliability_index_from_probability(1e-4)), 3)
    3.719
    """
    _interval("pf", pf, 0, 1, closed=False)
    return -special.ndtri(np.asarray(pf, dtype=float))
