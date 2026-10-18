r"""Open-channel flow: Manning's formula for uniform flow.

In a long channel of constant slope the flow settles at the normal depth,
where gravity along the bed balances friction on the wetted perimeter.
Chézy (1769) wrote :math:`V = C\sqrt{RS}`; Manning (1889) found that
:math:`C` grows like :math:`R^{1/6}` and gave the formula used today,
:math:`V = R^{2/3}S^{1/2}/n` in SI units.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["hydraulic_radius", "manning_velocity", "manning_flow"]


def hydraulic_radius(area, wetted_perimeter):
    r"""Hydraulic radius :math:`R = A/P` of a channel section.

    Parameters
    ----------
    area : float or array_like
        Flow area [m^2], positive.
    wetted_perimeter : float or array_like
        Wetted perimeter, excluding the free surface [m], positive.

    Returns
    -------
    float or ndarray
        Hydraulic radius [m]; a quarter of the hydraulic diameter.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    A rectangular channel 4 m wide flowing 1 m deep:

    >>> hydraulic_radius(4.0, 6.0)
    0.6666666666666666
    """
    _finite(area=area, wetted_perimeter=wetted_perimeter)
    _positive(area=area, wetted_perimeter=wetted_perimeter)
    return area / wetted_perimeter


def manning_velocity(n, R, slope):
    r"""Mean velocity of uniform open-channel flow, :math:`V = R^{2/3}S^{1/2}/n` (SI).

    Parameters
    ----------
    n : float or array_like
        Manning roughness coefficient [s/m^(1/3)], positive (about 0.013
        for finished concrete, 0.035 for a natural stream).
    R : float or array_like
        Hydraulic radius [m], positive.
    slope : float or array_like
        Bed slope :math:`S` [-], nonnegative (equal to the energy slope in
        uniform flow).

    Returns
    -------
    float or ndarray
        Mean velocity [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The coefficient 1 holds in SI only;
    US customary units need 1.486. Fully rough turbulent flow is assumed.

    References
    ----------
    R. Manning, "On the flow of water in open channels and pipes," *Trans.
    Inst. Civil Eng. Ireland* 20 (1891) 161-207 (read 4 December 1889).

    Examples
    --------
    >>> round(float(manning_velocity(0.013, 1.0, 0.001)), 3)
    2.433
    """
    _finite(n=n, R=R, slope=slope)
    _positive(n=n, R=R)
    _nonnegative(slope=slope)
    return np.asarray(R, dtype=float) ** (2 / 3) * np.sqrt(slope) / n


def manning_flow(n, area, wetted_perimeter, slope):
    r"""Volume flow of uniform open-channel flow, :math:`Q = A R^{2/3}S^{1/2}/n` (SI).

    Parameters
    ----------
    n : float or array_like
        Manning roughness coefficient [s/m^(1/3)], positive.
    area : float or array_like
        Flow area [m^2], positive.
    wetted_perimeter : float or array_like
        Wetted perimeter [m], positive.
    slope : float or array_like
        Bed slope [-], nonnegative.

    Returns
    -------
    float or ndarray
        Volume flow [m^3/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Solve this for the depth to get the
    normal depth of a channel.

    Examples
    --------
    >>> round(float(manning_flow(0.013, 4.0, 6.0, 0.001)), 3)
    7.425
    """
    return area * manning_velocity(n, hydraulic_radius(area, wetted_perimeter), slope)
