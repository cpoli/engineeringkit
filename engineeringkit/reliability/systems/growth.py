r"""Reliability growth: Duane's learning curve.

During development, each failure found and fixed makes the next one rarer.
Duane observed that the *cumulative* MTBF, total test time divided by total
failures, then grows as a power of test time, a straight line on log-log
paper.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["reliability_growth_duane"]


def reliability_growth_duane(t, K, alpha):
    r"""Duane's cumulative MTBF :math:`\theta_c(t) = K t^{\alpha}` after test time ``t``.

    Parameters
    ----------
    t : float or array_like
        Cumulative test time [s], positive.
    K : float
        Coefficient [s^(1-alpha)], positive: the cumulative MTBF at
        ``t = 1``.
    alpha : float
        Growth rate [-], nonnegative; Duane found 0.3-0.6 for vigorous
        programmes.

    Returns
    -------
    float or ndarray
        Cumulative MTBF [s]. The instantaneous (current) MTBF is
        :math:`\theta_c/(1-\alpha)`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    J. T. Duane, "Learning curve approach to reliability monitoring," *IEEE
    Trans. Aerospace* 2 (1964) 563-566.

    Examples
    --------
    >>> float(reliability_growth_duane(4.0, 2.0, 0.5))
    4.0
    """
    _finite(t=t, K=K, alpha=alpha)
    _positive(t=t, K=K)
    _nonnegative(alpha=alpha)
    return K * np.asarray(t) ** alpha
