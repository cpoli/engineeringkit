"""Decibel conversions shared by the propagation, noise and link-budget models."""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive

__all__ = ["db", "from_db", "dbm", "from_dbm"]


def db(x):
    r"""Power ratio in decibels, :math:`10\log_{10} x`.

    Parameters
    ----------
    x : float or array_like
        Linear power ratio [-], positive.

    Returns
    -------
    float or ndarray
        Ratio [dB].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(db(100.0))
    20.0
    """
    _positive(x=x)
    return 10 * np.log10(x)


def from_db(x_db):
    r"""Linear power ratio from decibels, :math:`10^{x/10}`.

    Parameters
    ----------
    x_db : float or array_like
        Ratio [dB].

    Returns
    -------
    float or ndarray
        Linear power ratio [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(from_db(30.0))
    1000.0
    """
    _finite(x_db=x_db)
    return 10 ** (np.asarray(x_db) / 10)


def dbm(power_w):
    r"""Power in decibels relative to one milliwatt.

    Parameters
    ----------
    power_w : float or array_like
        Power [W], positive.

    Returns
    -------
    float or ndarray
        Power [dBm].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(dbm(1.0))
    30.0
    """
    _positive(power_w=power_w)
    return 10 * np.log10(np.asarray(power_w) / 1e-3)


def from_dbm(x):
    r"""Power in watts from dBm, :math:`10^{-3}\cdot 10^{x/10}`.

    Parameters
    ----------
    x : float or array_like
        Power [dBm].

    Returns
    -------
    float or ndarray
        Power [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(from_dbm(30.0))
    1.0
    """
    _finite(x=x)
    return 1e-3 * 10 ** (np.asarray(x) / 10)
