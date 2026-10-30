r"""Flow measurement: Pitot's tube, Venturi's contraction, orifices and weirs.

Both instruments turn a speed into a pressure difference through Bernoulli's
equation. Pitot's bent tube (1732) faces the flow and reads the stagnation
rise :math:`\tfrac12\rho v^2`. Venturi's contraction (1797), made into a
water meter by Herschel (1887), speeds the flow up in a throat, and the
pressure drop there, with continuity, gives the volume flow.

Torricelli (1644) found that water leaves a small hole at the speed it
would gain falling from the free surface, :math:`\sqrt{2gh}`; the jet
contracts past the hole, so the flow is a discharge coefficient times the
ideal one. A weir is the same idea for an open channel: the head over a
sharp crest gives the flow, as :math:`H^{3/2}` for a rectangular notch
and as :math:`H^{5/2}` for James Thomson's triangular notch (1858).
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G

__all__ = [
    "pitot_velocity",
    "venturi_flow",
    "orifice_flow",
    "tank_drain_time",
    "rectangular_weir_flow",
    "rehbock_discharge_coefficient",
    "v_notch_weir_flow",
]


def pitot_velocity(delta_p, rho):
    r"""Flow speed from a Pitot-static reading, :math:`v = \sqrt{2\Delta p/\rho}`.

    Parameters
    ----------
    delta_p : float or array_like
        Stagnation minus static pressure [Pa], nonnegative.
    rho : float or array_like
        Fluid density [kg/m^3], positive.

    Returns
    -------
    float or ndarray
        Flow speed [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Incompressible flow; for gases the
    error is below 1% up to about Mach 0.2.

    References
    ----------
    H. Pitot, "Description d'une machine pour mesurer la vitesse des eaux
    courantes et le sillage des vaisseaux," *Hist. Acad. R. Sci.* (1732)
    363-376.

    Examples
    --------
    >>> float(pitot_velocity(2000.0, 1000.0))
    2.0
    """
    _finite(delta_p=delta_p, rho=rho)
    _nonnegative(delta_p=delta_p)
    _positive(rho=rho)
    return np.sqrt(2 * delta_p / rho)


def venturi_flow(delta_p, rho, A1, A2, discharge_coefficient=1.0):
    r"""Volume flow through a Venturi meter.

    .. math::

        Q = C_d\,A_2\sqrt{\frac{2\Delta p}{\rho\,(1 - (A_2/A_1)^2)}}

    Parameters
    ----------
    delta_p : float or array_like
        Inlet minus throat pressure [Pa], nonnegative.
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    A1 : float or array_like
        Inlet area [m^2], positive.
    A2 : float or array_like
        Throat area [m^2], positive and smaller than ``A1``.
    discharge_coefficient : float or array_like, optional
        :math:`C_d` in (0, 1]; the default 1 is the ideal meter, and 0.98 is
        typical of a machined Venturi tube.

    Returns
    -------
    float or ndarray
        Volume flow [m^3/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Horizontal meter and incompressible
    flow.

    References
    ----------
    G. B. Venturi, *Recherches expérimentales sur le principe de la
    communication latérale du mouvement dans les fluides* (Paris, 1797);
    C. Herschel, "The Venturi water meter," *Trans. ASCE* 17 (1887) 228-258.

    Examples
    --------
    >>> round(float(venturi_flow(15000.0, 1000.0, 0.02, 0.01)), 4)
    0.0632
    """
    _finite(delta_p=delta_p, rho=rho, A1=A1, A2=A2)
    _nonnegative(delta_p=delta_p)
    _positive(rho=rho, A1=A1, A2=A2)
    _interval("discharge_coefficient", discharge_coefficient, 0, 1, closed="right")
    if np.any(np.asarray(A2) >= np.asarray(A1)):
        raise ValueError("throat A2 must be smaller than A1")
    return discharge_coefficient * A1 * A2 * np.sqrt(2 * delta_p / (rho * (A1 * A1 - A2 * A2)))


def orifice_flow(area, head, discharge_coefficient=0.61, g=G):
    r"""Flow through a small orifice by Torricelli's law, :math:`Q = C_d A\sqrt{2gh}`.

    Parameters
    ----------
    area : float or array_like
        Orifice area [m^2], positive.
    head : float or array_like
        Height of the free surface above the orifice centre [m],
        nonnegative; a pressure difference :math:`\Delta p` across the
        orifice counts as :math:`\Delta p/(\rho g)`.
    discharge_coefficient : float or array_like, optional
        :math:`C_d` in (0, 1]; the default 0.61 is a sharp-edged orifice,
        whose jet contracts to about 61% of the hole area (the vena
        contracta), and 1 is the ideal Torricelli jet.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Volume flow [m^3/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The orifice is small next to the
    head and the approach velocity is negligible. An orifice plate in a
    pipe, where it is not, is :func:`venturi_flow` with its own
    :math:`C_d` (about 0.6).

    References
    ----------
    E. Torricelli, "De motu aquarum," in *De motu gravium naturaliter
    descendentium et proiectorum*, *Opera geometrica* (Florence: Masse and
    de Landis, 1644) 191-203.

    Examples
    --------
    >>> float(orifice_flow(0.01, 5.0, 1.0, g=10.0))
    0.1
    """
    _finite(area=area, head=head, discharge_coefficient=discharge_coefficient, g=g)
    _positive(area=area, g=g)
    _nonnegative(head=head)
    _interval("discharge_coefficient", discharge_coefficient, 0, 1, closed="right")
    return discharge_coefficient * area * np.sqrt(2 * g * head)


def tank_drain_time(tank_area, orifice_area, h0, h1=0.0, discharge_coefficient=0.61, g=G):
    r"""Time for a tank of constant section to drain from ``h0`` to ``h1`` through an orifice.

    .. math::

        t = \frac{A_t}{C_d A_o}\sqrt{\frac{2}{g}}\left(\sqrt{h_0} - \sqrt{h_1}\right)

    Parameters
    ----------
    tank_area : float or array_like
        Horizontal area of the tank [m^2], positive.
    orifice_area : float or array_like
        Orifice area [m^2], positive and much smaller than the tank.
    h0 : float or array_like
        Starting head above the orifice [m], nonnegative.
    h1 : float or array_like, optional
        Final head [m], in [0, ``h0``].
    discharge_coefficient : float or array_like, optional
        :math:`C_d` in (0, 1].
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Drain time [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Integrates :math:`A_t\,dh/dt =
    -Q(h)` from :func:`orifice_flow`, treating the flow as quasi-steady.
    The level falls fastest at the start, so the first half of the depth
    drains in less than a third of the time.

    Examples
    --------
    >>> round(float(tank_drain_time(1.0, 0.01, 5.0, discharge_coefficient=1.0, g=10.0)), 6)
    100.0
    """
    _finite(tank_area=tank_area, orifice_area=orifice_area, h0=h0, h1=h1, discharge_coefficient=discharge_coefficient, g=g)
    _positive(tank_area=tank_area, orifice_area=orifice_area, g=g)
    _nonnegative(h1=h1)
    _interval("discharge_coefficient", discharge_coefficient, 0, 1, closed="right")
    if np.any(np.asarray(h1) > np.asarray(h0)):
        raise ValueError("h1 must not exceed h0")
    return tank_area / (discharge_coefficient * orifice_area) * np.sqrt(2 / g) * (np.sqrt(h0) - np.sqrt(h1))


def rectangular_weir_flow(head, width, discharge_coefficient=0.62, g=G):
    r"""Flow over a sharp-crested rectangular weir, :math:`Q = C_d\,\tfrac23\,b\sqrt{2g}\,H^{3/2}`.

    Parameters
    ----------
    head : float or array_like
        Head over the crest, measured upstream of the drawdown [m],
        nonnegative.
    width : float or array_like
        Crest width [m], positive.
    discharge_coefficient : float or array_like, optional
        :math:`C_d` in (0, 1]; about 0.62, rising slowly with
        :math:`H/P` (see :func:`rehbock_discharge_coefficient`).
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Volume flow [m^3/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The ideal flow integrates
    Torricelli's speed :math:`\sqrt{2gz}` over the depth below the free
    surface. A full-width weir with a ventilated nappe; end contractions
    reduce the effective width by about :math:`0.1H` per end (Francis).

    Examples
    --------
    >>> round(float(rectangular_weir_flow(0.3, 2.0)), 4)
    0.6016
    """
    _finite(head=head, width=width, discharge_coefficient=discharge_coefficient, g=g)
    _nonnegative(head=head)
    _positive(width=width, g=g)
    _interval("discharge_coefficient", discharge_coefficient, 0, 1, closed="right")
    return discharge_coefficient * 2 / 3 * width * np.sqrt(2 * g) * np.asarray(head, dtype=float) ** 1.5


def rehbock_discharge_coefficient(head, weir_height):
    r"""Rehbock's discharge coefficient of a rectangular weir, :math:`C_d = 0.611 + 0.075\,H/P`.

    Parameters
    ----------
    head : float or array_like
        Head over the crest [m], nonnegative.
    weir_height : float or array_like
        Height of the crest above the channel bed :math:`P` [m], positive.

    Returns
    -------
    float or ndarray
        Discharge coefficient [-] for :func:`rectangular_weir_flow`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The :math:`H/P` term accounts for
    the approach velocity. Fitted for :math:`H/P` up to about 5, as given
    in Munson et al., *Fundamentals of Fluid Mechanics*.

    Examples
    --------
    >>> rehbock_discharge_coefficient(0.5, 1.0)
    0.6485
    """
    _finite(head=head, weir_height=weir_height)
    _nonnegative(head=head)
    _positive(weir_height=weir_height)
    return 0.611 + 0.075 * head / weir_height


def v_notch_weir_flow(head, angle=np.pi / 2, discharge_coefficient=0.58, g=G):
    r"""Flow over a sharp-crested triangular (V-notch) weir.

    .. math::

        Q = C_d\,\frac{8}{15}\,\sqrt{2g}\,\tan\frac{\theta}{2}\,H^{5/2}

    Parameters
    ----------
    head : float or array_like
        Head above the vertex of the notch [m], nonnegative.
    angle : float or array_like, optional
        Included angle of the notch :math:`\theta` [rad], in (0, π);
        the default is the 90° notch.
    discharge_coefficient : float or array_like, optional
        :math:`C_d` in (0, 1]; about 0.58 for a 90° notch.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Volume flow [m^3/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The flow section stays similar as
    the head changes, so :math:`C_d` is nearly constant, and the
    :math:`H^{5/2}` law keeps the head readable at small flows.

    References
    ----------
    J. Thomson, "On experiments on the measurement of water by triangular
    notches in weir boards," *Rep. Br. Assoc. Adv. Sci.* (1858), Trans.
    Sect., 181; (1861), Trans. Sect., 151.

    Examples
    --------
    >>> round(float(v_notch_weir_flow(0.2)), 5)
    0.02451
    """
    _finite(head=head, angle=angle, discharge_coefficient=discharge_coefficient, g=g)
    _nonnegative(head=head)
    _positive(g=g)
    _interval("angle", angle, 0, np.pi, closed=False)
    _interval("discharge_coefficient", discharge_coefficient, 0, 1, closed="right")
    return discharge_coefficient * 8 / 15 * np.sqrt(2 * g) * np.tan(angle / 2) * np.asarray(head, dtype=float) ** 2.5
