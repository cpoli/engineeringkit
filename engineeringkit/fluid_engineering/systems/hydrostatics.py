r"""Hydrostatics: Archimedes's buoyancy and Stevin's pressure on a submerged surface.

In a fluid at rest the pressure grows linearly with depth,
:math:`p = \rho g h`, whatever the shape of the vessel (Stevin, 1586). The
force on a plane surface is that pressure at the centroid times the area,
and it acts lower down, at the centre of pressure. A body immersed in the
fluid feels the net upward force :math:`\rho g V` (Archimedes).
"""

from __future__ import annotations

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G

__all__ = [
    "hydrostatic_pressure",
    "pressure_head",
    "buoyancy_force",
    "hydrostatic_force_plane",
    "center_of_pressure_depth",
]


def hydrostatic_pressure(rho, h, g=G):
    r"""Gauge pressure at depth :math:`h` in a fluid at rest, :math:`p = \rho g h`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    h : float or array_like
        Depth below the free surface [m], nonnegative.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Gauge pressure [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The density is constant
    (incompressible liquid).

    References
    ----------
    S. Stevin, *De Beghinselen des Waterwichts* (Leiden: Plantijn, 1586).

    Examples
    --------
    >>> round(hydrostatic_pressure(1000.0, 10.0), 1)
    98066.5
    """
    _finite(rho=rho, h=h, g=g)
    _positive(rho=rho, g=g)
    _nonnegative(h=h)
    return rho * g * h


def pressure_head(p, rho, g=G):
    r"""Height of a fluid column that produces pressure :math:`p`, :math:`h = p/(\rho g)`.

    Parameters
    ----------
    p : float or array_like
        Pressure [Pa]; gauge or absolute, and the head is of the same kind.
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Pressure head [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> round(pressure_head(101325.0, 13595.1), 4)
    0.76
    """
    _finite(p=p, rho=rho, g=g)
    _positive(rho=rho, g=g)
    return p / (rho * g)


def buoyancy_force(rho, displaced_volume, g=G):
    r"""Archimedes's buoyant force, :math:`F_B = \rho g V`.

    Parameters
    ----------
    rho : float or array_like
        Density of the surrounding fluid [kg/m^3], positive.
    displaced_volume : float or array_like
        Volume of fluid displaced by the body [m^3], nonnegative.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Upward force [N], equal to the weight of the displaced fluid.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    References
    ----------
    Archimedes, *On Floating Bodies*, Book I, Prop. 7 (c. 250 BC), in
    T. L. Heath, *The Works of Archimedes* (Cambridge, 1897).

    Examples
    --------
    >>> round(buoyancy_force(1000.0, 0.001), 3)
    9.807
    """
    _finite(rho=rho, displaced_volume=displaced_volume, g=g)
    _positive(rho=rho, g=g)
    _nonnegative(displaced_volume=displaced_volume)
    return rho * g * displaced_volume


def hydrostatic_force_plane(rho, depth_centroid, area, g=G):
    r"""Resultant pressure force on a submerged plane surface, :math:`F = \rho g h_c A`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    depth_centroid : float or array_like
        Vertical depth of the surface's centroid below the free surface [m],
        nonnegative.
    area : float or array_like
        Wetted area of the surface [m^2], positive.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Force normal to the surface [N], from gauge pressure only.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The result holds for any orientation
    of the plane, because the mean pressure over a plane is the pressure at
    its centroid.

    References
    ----------
    S. Stevin, *De Beghinselen des Waterwichts* (Leiden: Plantijn, 1586).

    Examples
    --------
    >>> round(hydrostatic_force_plane(1000.0, 1.0, 2.0))
    19613
    """
    _finite(rho=rho, depth_centroid=depth_centroid, area=area, g=g)
    _positive(rho=rho, area=area, g=g)
    _nonnegative(depth_centroid=depth_centroid)
    return rho * g * depth_centroid * area


def center_of_pressure_depth(centroid_distance, area, second_moment):
    r"""Distance to the centre of pressure on a submerged plane, :math:`y_p = y_c + I_c/(y_c A)`.

    Distances are measured along the plane from the line where its plane
    meets the free surface; for a vertical surface they are depths.

    Parameters
    ----------
    centroid_distance : float or array_like
        Distance :math:`y_c` of the centroid from the free-surface line [m], positive.
    area : float or array_like
        Wetted area [m^2], positive.
    second_moment : float or array_like
        Second moment of area about the horizontal centroidal axis [m^4],
        nonnegative.

    Returns
    -------
    float or ndarray
        Distance :math:`y_p` of the centre of pressure [m]; always below the centroid.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Gauge pressure only: the atmosphere
    acting on both faces adds no moment.

    Examples
    --------
    A vertical rectangle 1 m wide reaching 3 m down from the surface:

    >>> center_of_pressure_depth(1.5, 3.0, 3.0**3 / 12)
    2.0
    """
    _finite(centroid_distance=centroid_distance, area=area, second_moment=second_moment)
    _positive(centroid_distance=centroid_distance, area=area)
    _nonnegative(second_moment=second_moment)
    return centroid_distance + second_moment / (centroid_distance * area)
