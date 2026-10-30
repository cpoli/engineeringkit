r"""Section properties: second moment, elastic and plastic section moduli, built-up and transformed sections.

For a rectangle of width :math:`b` and depth :math:`h`, the second moment is
:math:`bh^3/12`, the elastic modulus :math:`Z = bh^2/6` (first yield at
:math:`M_y = \sigma_y Z`) and the plastic modulus :math:`Z_p = bh^2/4`
(full plastification at :math:`M_p = \sigma_y Z_p`). I, T, channel, angle
and hollow box sections are assembled from rectangles with the parallel-axis
theorem; weighting a rectangle by a modular ratio gives the transformed
section of a composite beam.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.structural_engineering.core.base import SectionProperties

__all__ = [
    "rectangular_second_moment",
    "rectangular_section_modulus",
    "plastic_section_modulus",
    "parallel_axis",
    "composite_section",
    "i_section",
    "t_section",
    "channel_section",
    "angle_section",
    "box_section",
]


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


def composite_section(rectangles):
    r"""Area, centroid and second moments of a section made of rectangles.

    Each part adds its own centroidal second moment and, by the parallel-axis
    theorem, its area times the squared distance to the section centroid:

    .. math::

        I_z = \sum n_i\left(\frac{b_ih_i^3}{12} + b_ih_i(y_i - \bar y)^2\right),\qquad
        I_{yz} = \sum n_ib_ih_i(y_i - \bar y)(z_i - \bar z).

    A weight :math:`n_i = -1` removes a hole; a modular ratio
    :math:`n_i = E_i/E_\mathrm{ref}` gives the transformed section of a
    composite beam, whose bending stress in part :math:`i` is
    :math:`n_i` times the flexure-formula stress of the transformed section.

    Parameters
    ----------
    rectangles : sequence of tuple
        ``(b, h, z, y)`` or ``(b, h, z, y, weight)`` for each part: width
        along :math:`z` [m] and height along :math:`y` [m], both positive,
        the coordinates of its centre [m], and its weight [-], nonzero
        (1 by default). Holes must lie inside solid parts.

    Returns
    -------
    SectionProperties

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    R. C. Hibbeler, *Mechanics of Materials*, 10th ed. (Pearson, 2017),
    section 6.6 and appendix A.
    """
    if not len(rectangles):
        raise ValueError("at least one rectangle is required")
    parts = []
    for part in rectangles:
        if len(part) not in (4, 5):
            raise ValueError("each rectangle must be (b, h, z, y) or (b, h, z, y, weight)")
        b, h, z, y = part[:4]
        weight = part[4] if len(part) == 5 else 1.0
        _positive(b=b, h=h)
        _finite(z=z, y=y, weight=weight)
        if weight == 0:
            raise ValueError("weight must be nonzero")
        parts.append((float(b), float(h), float(z), float(y), float(weight)))
    b, h, z, y, n = np.array(parts).T
    A = n * b * h
    area = float(A.sum())
    if area <= 0:
        raise ValueError("the net (transformed) area must be positive")
    yc, zc = float(A @ y / area), float(A @ z / area)
    Iz = float(np.sum(n * b * h**3 / 12 + A * (y - yc) ** 2))
    Iy = float(np.sum(n * h * b**3 / 12 + A * (z - zc) ** 2))
    Iyz = float(np.sum(A * (y - yc) * (z - zc)))
    return SectionProperties(area, yc, zc, Iz, Iy, Iyz, tuple(parts))


def i_section(depth, flange_width, flange_thickness, web_thickness):
    r"""Properties of a doubly symmetric I (wide-flange) section.

    Parameters
    ----------
    depth : float
        Overall depth [m], more than twice ``flange_thickness``.
    flange_width : float
        Flange width [m], at least ``web_thickness``.
    flange_thickness : float
        Flange thickness [m], positive.
    web_thickness : float
        Web thickness [m], positive.

    Returns
    -------
    SectionProperties
        Bottom of the section at :math:`y = 0`, web centred on :math:`z = 0`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Root fillets are ignored.
    """
    _positive(depth=depth, flange_width=flange_width, flange_thickness=flange_thickness, web_thickness=web_thickness)
    _positive(web_depth=depth - 2 * flange_thickness)
    _nonnegative(flange_outstand=flange_width - web_thickness)
    d, bf, tf, tw = depth, flange_width, flange_thickness, web_thickness
    return composite_section([(bf, tf, 0.0, tf / 2), (tw, d - 2 * tf, 0.0, d / 2), (bf, tf, 0.0, d - tf / 2)])


def t_section(depth, flange_width, flange_thickness, web_thickness):
    r"""Properties of a T section, flange on top.

    Parameters
    ----------
    depth : float
        Overall depth [m], more than ``flange_thickness``.
    flange_width : float
        Flange width [m], at least ``web_thickness``.
    flange_thickness : float
        Flange thickness [m], positive.
    web_thickness : float
        Web (stem) thickness [m], positive.

    Returns
    -------
    SectionProperties
        Bottom of the stem at :math:`y = 0`, stem centred on :math:`z = 0`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(depth=depth, flange_width=flange_width, flange_thickness=flange_thickness, web_thickness=web_thickness)
    _positive(stem_depth=depth - flange_thickness)
    _nonnegative(flange_outstand=flange_width - web_thickness)
    d, bf, tf, tw = depth, flange_width, flange_thickness, web_thickness
    return composite_section([(tw, d - tf, 0.0, (d - tf) / 2), (bf, tf, 0.0, d - tf / 2)])


def channel_section(depth, flange_width, flange_thickness, web_thickness):
    r"""Properties of a channel (C) section, flanges pointing in :math:`+z`.

    Parameters
    ----------
    depth : float
        Overall depth [m], more than twice ``flange_thickness``.
    flange_width : float
        Overall flange width, including the web [m], more than ``web_thickness``.
    flange_thickness : float
        Flange thickness [m], positive.
    web_thickness : float
        Web thickness [m], positive.

    Returns
    -------
    SectionProperties
        Back of the web at :math:`z = 0`, bottom at :math:`y = 0`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. A channel loaded through its centroid twists:
    see :func:`channel_shear_center`.
    """
    _positive(depth=depth, flange_width=flange_width, flange_thickness=flange_thickness, web_thickness=web_thickness)
    _positive(web_depth=depth - 2 * flange_thickness, flange_outstand=flange_width - web_thickness)
    d, bf, tf, tw = depth, flange_width, flange_thickness, web_thickness
    zf = tw + (bf - tw) / 2
    return composite_section([(tw, d, tw / 2, d / 2), (bf - tw, tf, zf, tf / 2), (bf - tw, tf, zf, d - tf / 2)])


def angle_section(leg_height, leg_width, thickness):
    r"""Properties of an angle (L) section, the simplest section with :math:`I_{yz} \ne 0`.

    Parameters
    ----------
    leg_height : float
        Overall height of the vertical leg [m], more than ``thickness``.
    leg_width : float
        Overall width of the horizontal leg [m], more than ``thickness``.
    thickness : float
        Leg thickness [m], positive.

    Returns
    -------
    SectionProperties
        Heel (outer corner) at the origin, legs along :math:`+y` and :math:`+z`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(leg_height=leg_height, leg_width=leg_width, thickness=thickness)
    _positive(vertical_outstand=leg_height - thickness, horizontal_outstand=leg_width - thickness)
    h, b, t = leg_height, leg_width, thickness
    return composite_section([(t, h, t / 2, h / 2), (b - t, t, t + (b - t) / 2, t / 2)])


def box_section(depth, width, thickness):
    r"""Properties of a rectangular hollow (box) section of uniform wall thickness.

    Parameters
    ----------
    depth : float
        Overall depth [m], more than twice ``thickness``.
    width : float
        Overall width [m], more than twice ``thickness``.
    thickness : float
        Wall thickness [m], positive.

    Returns
    -------
    SectionProperties
        Bottom at :math:`y = 0`, centred on :math:`z = 0`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Corner radii are ignored.
    """
    _positive(depth=depth, width=width, thickness=thickness)
    _positive(inner_depth=depth - 2 * thickness, inner_width=width - 2 * thickness)
    d, b, t = depth, width, thickness
    return composite_section([(b, d, 0.0, d / 2), (b - 2 * t, d - 2 * t, 0.0, d / 2, -1.0)])
