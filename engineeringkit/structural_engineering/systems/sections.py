r"""Section properties: second moment, elastic and plastic section moduli.

For a rectangle of width :math:`b` and depth :math:`h`, the second moment is
:math:`bh^3/12`, the elastic modulus :math:`Z = bh^2/6` (first yield at
:math:`M_y = \sigma_y Z`) and the plastic modulus :math:`Z_p = bh^2/4`
(full plastification at :math:`M_p = \sigma_y Z_p`).
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["rectangular_second_moment", "rectangular_section_modulus", "plastic_section_modulus", "parallel_axis"]


def rectangular_second_moment(b, h):
    r"""Second moment of area of a rectangle about its centroidal axis, :math:`bh^3/12`.

    Parameters
    ----------
    b : float or array_like
        Width, parallel to the bending axis [m], positive.
    h : float or array_like
        Depth, in the plane of bending [m], positive.

    Returns
    -------
    float or ndarray
        Second moment of area [m^4].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(b=b, h=h)
    _positive(b=b, h=h)
    return np.asarray(b) * np.asarray(h) ** 3 / 12


def rectangular_section_modulus(b, h):
    r"""Elastic section modulus of a rectangle, :math:`Z = I/c = bh^2/6`.

    Parameters
    ----------
    b : float or array_like
        Width [m], positive.
    h : float or array_like
        Depth [m], positive.

    Returns
    -------
    float or ndarray
        Section modulus [m^3]; the extreme-fibre stress is :math:`M/Z`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    G. Galilei, *Discorsi e dimostrazioni matematiche intorno a due nuove
    scienze* (Leiden: Elsevier, 1638), second day.
    """
    _finite(b=b, h=h)
    _positive(b=b, h=h)
    return np.asarray(b) * np.asarray(h) ** 2 / 6


def plastic_section_modulus(b, h):
    r"""Plastic section modulus of a rectangle, :math:`Z_p = bh^2/4`.

    When the whole section has yielded, half of it is at :math:`+\sigma_y`
    and half at :math:`-\sigma_y`, and the moment is
    :math:`M_p = \sigma_y Z_p`. The ratio :math:`Z_p/Z = 1.5` is the
    rectangle's *shape factor*.

    Parameters
    ----------
    b : float or array_like
        Width [m], positive.
    h : float or array_like
        Depth [m], positive.

    Returns
    -------
    float or ndarray
        Plastic section modulus [m^3].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(b=b, h=h)
    _positive(b=b, h=h)
    return np.asarray(b) * np.asarray(h) ** 2 / 4


def parallel_axis(Ic, A, d):
    r"""Second moment about a parallel axis, :math:`I = I_c + Ad^2`.

    Parameters
    ----------
    Ic : float or array_like
        Second moment about the centroidal axis [m^4], nonnegative.
    A : float or array_like
        Area [m^2], positive.
    d : float or array_like
        Distance between the axes [m].

    Returns
    -------
    float or ndarray
        Second moment about the shifted axis [m^4].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Ic=Ic, A=A, d=d)
    _nonnegative(Ic=Ic)
    _positive(A=A)
    return Ic + np.asarray(A) * np.asarray(d) ** 2
