"""Calorically-perfect-gas helpers shared by the compressible-flow models."""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative

__all__ = ["check_gamma", "isentropic_factor"]


def check_gamma(gamma):
    """Require a ratio of specific heats ``gamma > 1`` (scalar or array).

    Parameters
    ----------
    gamma : float or array_like
        Ratio of specific heats :math:`c_p/c_v` [-].

    Raises
    ------
    ValueError
        If any value is nonfinite or not greater than one.
    """
    _interval("gamma", gamma, 1, np.inf, closed=False)


def isentropic_factor(M, gamma=1.4):
    r"""Stagnation-to-static temperature ratio :math:`1 + \tfrac{\gamma-1}{2}M^2`.

    Every isentropic relation (stagnation temperature and pressure, the
    area-Mach relation) is a power of this one factor.

    Parameters
    ----------
    M : float or array_like
        Mach number [-], nonnegative.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        :math:`T_0/T` [-].

    Examples
    --------
    >>> float(isentropic_factor(1.0))
    1.2
    """
    _nonnegative(M=M)
    check_gamma(gamma)
    return 1 + (gamma - 1) / 2 * np.asarray(M) ** 2
