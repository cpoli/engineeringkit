r"""Open-channel flow: Manning's formula, specific energy and the hydraulic jump.

In a long channel of constant slope the flow settles at the normal depth,
where gravity along the bed balances friction on the wetted perimeter.
Chézy (1769) wrote :math:`V = C\sqrt{RS}`; Manning (1889) found that
:math:`C` grows like :math:`R^{1/6}` and gave the formula used today,
:math:`V = R^{2/3}S^{1/2}/n` in SI units.

Over a short reach the bed slope and friction drop out and the energy
above the bed, the specific energy :math:`E = y + q^2/(2gy^2)` of
Bakhmeteff (1912), sets the depth: it is least at the critical depth, and
every larger energy can be carried at two alternate depths. A fast,
shallow (supercritical) stream that must slow down does so in a hydraulic
jump, whose downstream depth Bélanger (1828) found from the momentum
balance.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G

__all__ = [
    "hydraulic_radius",
    "manning_velocity",
    "manning_flow",
    "specific_energy",
    "critical_depth",
    "alternate_depths",
    "sequent_depth",
    "hydraulic_jump_head_loss",
]


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


def specific_energy(y, q, g=G):
    r"""Specific energy of flow in a wide rectangular channel, :math:`E = y + q^2/(2gy^2)`.

    Parameters
    ----------
    y : float or array_like
        Flow depth [m], positive.
    q : float or array_like
        Flow per unit width :math:`Q/b` [m^2/s], nonnegative.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Head above the channel bed [m]: the depth plus the velocity head.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Hydrostatic pressure (straight,
    parallel streamlines) and a uniform velocity over the depth.

    References
    ----------
    B. A. Bakhmeteff, *O neravnomernom dvizhenii zhidkosti v otkrytom
    rusle* [Varied flow of liquids in open channels] (St Petersburg, 1912);
    B. A. Bakhmeteff, *Hydraulics of Open Channels* (New York: McGraw-Hill,
    1932).

    Examples
    --------
    >>> specific_energy(1.0, 2.0, g=10.0)
    1.2
    """
    _finite(y=y, q=q, g=g)
    _positive(y=y, g=g)
    _nonnegative(q=q)
    return y + q * q / (2 * g * y * y)


def critical_depth(q, g=G):
    r"""Critical depth of a wide rectangular channel, :math:`y_c = (q^2/g)^{1/3}`.

    Parameters
    ----------
    q : float or array_like
        Flow per unit width [m^2/s], nonnegative.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Critical depth [m]: the depth of least specific energy,
        :math:`E_{min} = 3y_c/2`, where the Froude number is 1.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> round(float(critical_depth(np.sqrt(10.0), g=10.0)), 6)
    1.0
    """
    _finite(q=q, g=g)
    _nonnegative(q=q)
    _positive(g=g)
    return np.cbrt(q * q / g)


def alternate_depths(E, q, g=G):
    r"""The subcritical and supercritical depths that carry ``q`` at specific energy ``E``.

    Parameters
    ----------
    E : float
        Specific energy [m], at least :math:`3y_c/2`.
    q : float
        Flow per unit width [m^2/s], positive.
    g : float, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    tuple of float
        ``(y_subcritical, y_supercritical)`` [m]; both equal the critical
        depth when ``E`` is the minimum.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Scalar inputs only; each root is
    bracketed on its own side of the critical depth and found by Brent's
    method. A sluice gate or a step in the bed moves the flow between the
    two depths at constant :math:`E`.

    Examples
    --------
    >>> [round(y, 4) for y in alternate_depths(1.2, 2.0, g=10.0)]
    [1.0, 0.5583]
    """
    _finite(E=E, q=q, g=g)
    _positive(q=q, g=g)
    yc = float(critical_depth(q, g))
    if E < 1.5 * yc * (1 - 1e-12):
        raise ValueError(f"E must be at least the minimum specific energy {1.5 * yc:g} m")
    if E <= 1.5 * yc:
        return yc, yc

    def residual(y):
        return y + q * q / (2 * g * y * y) - E

    deep = optimize.brentq(residual, yc, E, xtol=1e-14)
    shallow = optimize.brentq(residual, q / np.sqrt(2 * g * E), yc, xtol=1e-14)
    return float(deep), float(shallow)


def sequent_depth(y1, v1, g=G):
    r"""Depth after a hydraulic jump from Bélanger's momentum balance.

    .. math::

        \frac{y_2}{y_1} = \frac{\sqrt{1 + 8\,\mathrm{Fr}_1^2} - 1}{2},
        \qquad \mathrm{Fr}_1 = \frac{v_1}{\sqrt{g y_1}}

    Parameters
    ----------
    y1 : float or array_like
        Depth before the jump [m], positive.
    v1 : float or array_like
        Mean velocity before the jump [m/s].
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Sequent (conjugate) depth [m]. A jump forms only from supercritical
        flow, :math:`\mathrm{Fr}_1 > 1`; the relation is symmetric, so a
        subcritical ``y1`` returns the supercritical depth paired with it.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Horizontal, wide rectangular channel
    with bed friction over the short jump neglected: the pressure force
    plus momentum flux :math:`\rho g y^2/2 + \rho q^2/y` per unit width is
    the same on both sides.

    References
    ----------
    J.-B. Bélanger, *Essai sur la solution numérique de quelques problèmes
    relatifs au mouvement permanent des eaux courantes* (Paris:
    Carilian-Goeury, 1828); H. Chanson, "Development of the Bélanger
    equation and backwater equation by Jean-Baptiste Bélanger (1828),"
    *J. Hydraul. Eng.* 135 (2009) 159-163,
    doi:10.1061/(ASCE)0733-9429(2009)135:3(159).

    Examples
    --------
    >>> round(float(sequent_depth(0.5, np.sqrt(15.0), g=10.0)), 6)
    1.0
    """
    _finite(y1=y1, v1=v1, g=g)
    _positive(y1=y1, g=g)
    Fr1 = v1 / np.sqrt(g * y1)
    return 0.5 * y1 * (np.sqrt(1 + 8 * Fr1 * Fr1) - 1)


def hydraulic_jump_head_loss(y1, y2):
    r"""Head lost in a hydraulic jump, :math:`h_L = (y_2 - y_1)^3/(4y_1y_2)`.

    Parameters
    ----------
    y1 : float or array_like
        Depth before the jump [m], positive.
    y2 : float or array_like
        Sequent depth after the jump [m], positive.

    Returns
    -------
    float or ndarray
        Head loss [m], the drop in specific energy across the jump.
        Negative (impossible) if ``y2 < y1``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Holds only for a pair of sequent
    depths from :func:`sequent_depth`.

    Examples
    --------
    >>> hydraulic_jump_head_loss(1.0, 3.0)
    0.6666666666666666
    """
    _finite(y1=y1, y2=y2)
    _positive(y1=y1, y2=y2)
    return (y2 - y1) ** 3 / (4 * y1 * y2)
