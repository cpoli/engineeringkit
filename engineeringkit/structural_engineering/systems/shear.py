r"""Transverse shear in beams: Jourawski's formula, shear flow and the shear centre.

Where the bending moment changes along a beam, the difference in normal
force on the part of the section beyond a cut must be carried by shear on
the cut. Per unit length this *shear flow* is :math:`q = VQ/I`, and the
average shear stress across a cut of width :math:`t` is
:math:`\tau = VQ/(It)`. In a thin-walled open section the flange shear flows
form a couple, so the resultant shear passes through the *shear centre*,
not the centroid: a channel loaded through its web twists.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["transverse_shear_stress", "shear_flow", "channel_shear_center"]


def transverse_shear_stress(V, Q, I, t):
    r"""Jourawski's average shear stress across a cut, :math:`\tau = VQ/(It)`.

    Parameters
    ----------
    V : float or array_like
        Shear force on the section [N].
    Q : float or array_like
        First moment, about the neutral axis, of the area beyond the cut
        [m^3] (:meth:`SectionProperties.first_moment`).
    I : float or array_like
        Second moment of the whole section about the neutral axis [m^4], positive.
    t : float or array_like
        Width of the cut [m], positive (:meth:`SectionProperties.width_at`).

    Returns
    -------
    float or ndarray
        Shear stress [Pa], averaged across the width of the cut. For a
        rectangle it peaks at :math:`1.5V/A` on the neutral axis.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. The stress is assumed uniform across the cut,
    which is accurate for narrow cuts and thin walls, not across the wide
    flange of a stocky section.

    References
    ----------
    D. J. Jourawski, "Sur la résistance d'un corps prismatique et d'une
    pièce composée en bois ou en tôle de fer à une force perpendiculaire à
    leur longueur," *Annales des ponts et chaussées*, 3e série, 12 (1856)
    328-351.
    """
    _finite(V=V, Q=Q)
    _positive(I=I, t=t)
    return np.asarray(V) * Q / (np.asarray(I) * t)


def shear_flow(V, Q, I):
    r"""Shear force per unit length on a longitudinal cut, :math:`q = VQ/I`.

    In a built-up beam this is the force the glue, nails or bolts joining
    the part beyond the cut must carry per metre; fasteners of capacity
    :math:`F` need a spacing of at most :math:`F/q`.

    Parameters
    ----------
    V : float or array_like
        Shear force on the section [N].
    Q : float or array_like
        First moment of the area beyond the cut about the neutral axis [m^3].
    I : float or array_like
        Second moment of the whole section [m^4], positive.

    Returns
    -------
    float or ndarray
        Shear flow [N/m].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(V=V, Q=Q)
    _positive(I=I)
    return np.asarray(V) * Q / I


def channel_shear_center(depth, flange_width, flange_thickness, web_thickness):
    r"""Distance of a thin-walled channel's shear centre behind its web.

    The flange shear flows, each with resultant :math:`F = Vb^2ht_f/(4I)`,
    form a couple :math:`Fh` that the shear force can only balance if it acts
    at

    .. math::

        e = \frac{3b^2t_f}{6bt_f + ht_w}

    from the web centreline, on the side away from the flanges.

    Parameters
    ----------
    depth : float
        Distance :math:`h` between the flange centrelines [m], positive.
    flange_width : float
        Distance :math:`b` from the web centreline to the flange tips [m],
        nonnegative.
    flange_thickness : float
        Flange thickness :math:`t_f` [m], positive.
    web_thickness : float
        Web thickness :math:`t_w` [m], positive.

    Returns
    -------
    float or ndarray
        Shear-centre offset :math:`e` [m]; :math:`3b^2/(6b + h)` for a
        uniform wall thickness.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Centreline (thin-wall) dimensions are used, and
    the web's contribution to the second moment is :math:`t_wh^3/12`.

    References
    ----------
    R. Maillart, "Zur Frage der Biegung," *Schweizerische Bauzeitung* 77
    (1921) 195-197.
    """
    _positive(depth=depth, flange_thickness=flange_thickness, web_thickness=web_thickness)
    _nonnegative(flange_width=flange_width)
    b, tf = np.asarray(flange_width, float), np.asarray(flange_thickness, float)
    return 3 * b**2 * tf / (6 * b * tf + np.asarray(depth) * web_thickness)
