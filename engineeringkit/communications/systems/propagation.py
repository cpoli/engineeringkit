r"""Free-space radio propagation and the decibel link budget.

Hertz's experiments made the wavelength :math:`\lambda = c/f` a measurable
quantity; Friis's transmission formula turned antenna gains, distance and
wavelength into received power; and the link budget adds the same terms in
decibels.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.communications.utils.decibels import db

__all__ = ["SPEED_OF_LIGHT", "wavelength", "free_space_path_loss", "free_space_path_loss_db", "friis_received_power", "link_budget_db"]

SPEED_OF_LIGHT = 299792458.0
"""Speed of light in vacuum [m/s] (exact)."""


def wavelength(f, c=SPEED_OF_LIGHT):
    r"""Free-space wavelength :math:`\lambda = c/f`.

    Parameters
    ----------
    f : float or array_like
        Frequency [Hz], positive.
    c : float, optional
        Propagation speed [m/s], positive.

    Returns
    -------
    float or ndarray
        Wavelength [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(wavelength(100e6)), 4)
    2.9979
    """
    _positive(f=f, c=c)
    return c / np.asarray(f)


def free_space_path_loss(distance, f, c=SPEED_OF_LIGHT):
    r"""Free-space path loss between isotropic antennas, :math:`(4\pi d f / c)^2`.

    Parameters
    ----------
    distance : float or array_like
        Separation [m], positive; far field is assumed.
    f : float or array_like
        Frequency [Hz], positive.
    c : float, optional
        Propagation speed [m/s].

    Returns
    -------
    float or ndarray
        Linear loss ratio [-], at least one in the far field.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(free_space_path_loss_db(1000.0, 2.4e9)), 2)
    100.05
    """
    _positive(distance=distance, f=f, c=c)
    return (4 * np.pi * np.asarray(distance) * f / c) ** 2


def free_space_path_loss_db(distance, f, c=SPEED_OF_LIGHT):
    """Free-space path loss in decibels; see :func:`free_space_path_loss`.

    Parameters
    ----------
    distance : float or array_like
        Separation [m], positive.
    f : float or array_like
        Frequency [Hz], positive.
    c : float, optional
        Propagation speed [m/s].

    Returns
    -------
    float or ndarray
        Path loss [dB].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(free_space_path_loss_db(1.0, 299792458.0 / (4 * np.pi))), 6)
    0.0
    """
    return db(free_space_path_loss(distance, f, c))


def friis_received_power(Pt, Gt, Gr, distance, f):
    r"""Friis transmission equation, :math:`P_r = P_t G_t G_r (\lambda/4\pi d)^2`.

    Parameters
    ----------
    Pt : float or array_like
        Transmitted power [W], nonnegative.
    Gt, Gr : float or array_like
        Transmit and receive antenna gains [-] (linear), positive.
    distance : float or array_like
        Separation [m], positive.
    f : float or array_like
        Frequency [Hz], positive.

    Returns
    -------
    float or ndarray
        Received power [W].

    Notes
    -----
    Free space, far-field antennas and matched polarization and impedance
    are assumed. Inputs are raw SI numbers; see :doc:`/guides/communications`.

    References
    ----------
    H. T. Friis, "A note on a simple transmission formula," *Proc. IRE* 34
    (1946) 254-256.

    Examples
    --------
    >>> round(float(friis_received_power(1.0, 1.0, 1.0, 1000.0, 2.4e9)) * 1e10, 4)
    0.9881
    """
    _nonnegative(Pt=Pt)
    _positive(Gt=Gt, Gr=Gr)
    return np.asarray(Pt) * Gt * Gr / free_space_path_loss(distance, f)


def link_budget_db(pt_dbm, gt_db, gr_db, path_loss_db, losses_db=0.0):
    r"""Received power from a decibel link budget, :math:`P_r = P_t + G_t + G_r - L_p - L`.

    Parameters
    ----------
    pt_dbm : float or array_like
        Transmitted power [dBm].
    gt_db, gr_db : float or array_like
        Antenna gains [dBi].
    path_loss_db : float or array_like
        Propagation loss [dB].
    losses_db : float or array_like, optional
        Other losses (cables, pointing, atmosphere) [dB].

    Returns
    -------
    float or ndarray
        Received power [dBm].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(link_budget_db(30.0, 10.0, 10.0, 100.0, 2.0))
    -52.0
    """
    _finite(pt_dbm=pt_dbm, gt_db=gt_db, gr_db=gr_db, path_loss_db=path_loss_db, losses_db=losses_db)
    return np.asarray(pt_dbm) + gt_db + gr_db - path_loss_db - losses_db
