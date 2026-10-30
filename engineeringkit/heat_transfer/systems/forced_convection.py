r"""Forced convection: flat plates, cylinders, tube banks, tubes and the heat and mass transfer analogy.

Graetz (1883) solved laminar heat transfer in a tube, the first convection
problem with an exact answer. Pohlhausen (1921) solved the thermal boundary
layer on a flat plate, and Colburn (1933) with Chilton (1934) carried heat
transfer results over to mass transfer and friction. Žukauskas (1972)
correlated tube banks in crossflow, and Churchill and Bernstein (1977) the
single cylinder over its whole Reynolds range.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "flat_plate_nusselt",
    "churchill_bernstein_nusselt",
    "tube_bank_max_velocity",
    "zukauskas_nusselt",
    "laminar_tube_nusselt",
    "entry_length",
    "hausen_nusselt",
    "lewis_number",
    "chilton_colburn_mass_transfer",
]


def flat_plate_nusselt(Re, Pr, Re_c=5e5, local=False):
    r"""Nusselt number of a flat plate in parallel flow, laminar, turbulent or mixed.

    Laminar (Pohlhausen): :math:`Nu_x = 0.332\,Re_x^{1/2}Pr^{1/3}` locally and
    :math:`\bar{Nu}_L = 0.664\,Re_L^{1/2}Pr^{1/3}` on average. Turbulent:
    :math:`Nu_x = 0.0296\,Re_x^{4/5}Pr^{1/3}`. Transition at ``Re_c`` makes
    the average over a plate with both regimes

    .. math::

        \bar{Nu}_L = \left[0.037\,(Re_L^{4/5} - Re_c^{4/5}) + 0.664\,Re_c^{1/2}\right] Pr^{1/3}.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number :math:`u_\infty x/\nu` on the distance from the leading
        edge (the plate length for the average) [-], nonnegative.
    Pr : float or array_like
        Prandtl number [-], positive; the :math:`Pr^{1/3}` laws are meant for
        :math:`0.6 \lesssim Pr \lesssim 60`.
    Re_c : float, optional
        Transition Reynolds number [-], nonnegative; 0 makes the boundary layer
        turbulent from the leading edge and ``numpy.inf`` keeps it laminar.
    local : bool, optional
        Return the local :math:`Nu_x = h_x x/k` instead of the average :math:`\bar h L/k`.

    Returns
    -------
    float or ndarray
        Nusselt number [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Isothermal plate, properties at the film
    temperature. A uniform heat flux raises the laminar coefficient 0.332
    to 0.453.

    References
    ----------
    E. Pohlhausen, "Der Wärmeaustausch zwischen festen Körpern und
    Flüssigkeiten mit kleiner Reibung und kleiner Wärmeleitung," *Z. Angew.
    Math. Mech.* 1 (1921) 115-121.

    Examples
    --------
    >>> round(float(flat_plate_nusselt(1e5, 0.7)), 1)
    186.4
    >>> round(float(flat_plate_nusselt(1e6, 0.7)), 1)
    1299.2
    """
    _nonnegative(Re=Re)
    _positive(Pr=Pr)
    if isinstance(Re_c, bool) or not np.ndim(Re_c) == 0 or not float(Re_c) >= 0:
        raise ValueError("Re_c must be a nonnegative number (inf allowed)")
    Re = np.asarray(Re, dtype=float)
    pr3 = np.asarray(Pr, dtype=float) ** (1 / 3)
    laminar = Re <= Re_c
    if local:
        return np.where(laminar, 0.332 * Re**0.5, 0.0296 * Re**0.8) * pr3
    Re_t = np.minimum(Re_c, Re)
    return np.where(laminar, 0.664 * Re**0.5, 0.037 * (Re**0.8 - Re_t**0.8) + 0.664 * Re_t**0.5) * pr3


def churchill_bernstein_nusselt(Re, Pr):
    r"""Churchill-Bernstein average Nusselt number for a circular cylinder in crossflow.

    .. math::

        \bar{Nu}_D = 0.3 + \frac{0.62\,Re_D^{1/2}Pr^{1/3}}{[1 + (0.4/Pr)^{2/3}]^{1/4}}
        \left[1 + \left(\frac{Re_D}{282\,000}\right)^{5/8}\right]^{4/5}

    Parameters
    ----------
    Re : float or array_like
        Reynolds number :math:`VD/\nu` on the diameter [-], nonnegative.
    Pr : float or array_like
        Prandtl number [-], positive.

    Returns
    -------
    float or ndarray
        Average Nusselt number :math:`\bar hD/k` [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Meant for :math:`Re\,Pr \gtrsim 0.2`;
    properties at the film temperature. One blended expression covers the
    whole range, from creeping flow to the drag crisis; the authors present it
    as a slight underestimate of most data.

    References
    ----------
    S. W. Churchill and M. Bernstein, "A correlating equation for forced
    convection from gases and liquids to a circular cylinder in crossflow,"
    *J. Heat Transfer* 99 (1977) 300-306.

    Examples
    --------
    >>> round(float(churchill_bernstein_nusselt(6071.0, 0.7)), 1)
    40.6
    """
    _nonnegative(Re=Re)
    _positive(Pr=Pr)
    Re, Pr = np.asarray(Re, dtype=float), np.asarray(Pr, dtype=float)
    laminar = 0.62 * Re**0.5 * Pr ** (1 / 3) / (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25
    return 0.3 + laminar * (1 + (Re / 282000) ** 0.625) ** 0.8


_ARRANGEMENTS = ("aligned", "staggered")


def tube_bank_max_velocity(V, D, ST, SL, arrangement="aligned"):
    r"""Maximum velocity between the tubes of a bank, on which Žukauskas's Reynolds number is based.

    Aligned: :math:`V_{max} = S_T V/(S_T - D)`. Staggered: the same, unless
    the diagonal gap :math:`S_D - D` with :math:`S_D = [S_L^2 + (S_T/2)^2]^{1/2}`
    is narrower than :math:`(S_T - D)/2`, in which case
    :math:`V_{max} = S_T V/[2(S_D - D)]`.

    Parameters
    ----------
    V : float or array_like
        Upstream (approach) velocity [m/s], nonnegative.
    D : float
        Tube outer diameter [m], positive.
    ST : float
        Transverse pitch, normal to the flow [m], greater than ``D``.
    SL : float
        Longitudinal pitch, along the flow [m], positive.
    arrangement : {"aligned", "staggered"}, optional
        Tube layout.

    Returns
    -------
    float or ndarray
        Maximum velocity [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(tube_bank_max_velocity(6.0, 0.01, 0.015, 0.015)), 6)
    18.0
    """
    _nonnegative(V=V)
    _positive(D=D, ST=ST, SL=SL)
    if arrangement not in _ARRANGEMENTS:
        raise ValueError(f"arrangement must be one of {list(_ARRANGEMENTS)}")
    if ST <= D:
        raise ValueError("ST must exceed D")
    gap = ST - D
    if arrangement == "staggered":
        SD = np.hypot(SL, ST / 2)
        if SD <= D:
            raise ValueError("diagonal pitch must exceed D")
        gap = min(gap, 2 * (SD - D))
    return ST * np.asarray(V, dtype=float) / gap


# Žukauskas (1972) as tabulated by Incropera et al., Table 7.5: (upper Re, C, m) per range.
_ZUKAUSKAS = {
    "aligned": ((1e2, 0.80, 0.40), (1e3, 0.52, 0.50), (2e5, 0.27, 0.63), (2e6, 0.021, 0.84)),
    "staggered": ((1e2, 0.90, 0.40), (1e3, 0.71, 0.50), (2e5, 0.35, 0.60), (2e6, 0.022, 0.84)),
}
# Row correction C2 for fewer than 20 rows (Incropera et al., Table 7.6), 1 at 20 rows.
_ROWS = (1, 2, 3, 4, 5, 7, 10, 13, 16, 20)
_ROW_FACTOR = {
    "aligned": (0.70, 0.80, 0.86, 0.90, 0.92, 0.95, 0.97, 0.98, 0.99, 1.0),
    "staggered": (0.64, 0.76, 0.84, 0.89, 0.92, 0.95, 0.97, 0.98, 0.99, 1.0),
}


def zukauskas_nusselt(Re, Pr, ST, SL, arrangement="aligned", rows=20, Pr_s=None):
    r"""Žukauskas's average Nusselt number for a bank of tubes in crossflow.

    .. math::

        \bar{Nu}_D = C_2\,C\,Re_{D,max}^{m}\,Pr^{0.36}\,(Pr/Pr_s)^{1/4},

    with :math:`C` and :math:`m` from Žukauskas's table for the layout and
    Reynolds range, and :math:`C_2` the correction for fewer than 20 rows.
    In the staggered range :math:`10^3 < Re < 2\times10^5`, :math:`C =
    0.35\,(S_T/S_L)^{1/5}` for :math:`S_T/S_L < 2` and 0.40 above.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number :math:`V_{max}D/\nu` on the maximum velocity (see
        :func:`tube_bank_max_velocity`) [-], in [10, 2e6].
    Pr : float or array_like
        Prandtl number at the mean fluid temperature [-], positive; meant for
        :math:`0.7 \lesssim Pr \lesssim 500`.
    ST, SL : float
        Transverse and longitudinal pitches [m], positive.
    arrangement : {"aligned", "staggered"}, optional
        Tube layout.
    rows : int, optional
        Number of tube rows along the flow, at least one; 20 or more needs no correction.
    Pr_s : float or array_like, optional
        Prandtl number at the tube surface temperature [-], positive; ``None``
        takes it equal to ``Pr``.

    Returns
    -------
    float or ndarray
        Average Nusselt number :math:`\bar hD/k` over the bank [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Properties except :math:`Pr_s` at the
    arithmetic mean of the inlet and outlet fluid temperatures. Between 100
    and 1000 the table values 0.52 and 0.71 with :math:`m = 0.5` are
    Žukauskas's; Incropera et al. suggest the single-cylinder correlation
    there. Aligned banks with :math:`S_T/S_L < 0.7` transfer heat poorly and
    are not covered by the correlation. The row correction is interpolated
    linearly between the tabulated row counts.

    References
    ----------
    A. Žukauskas, "Heat transfer from tubes in crossflow," *Adv. Heat
    Transfer* 8 (1972) 93-160.

    Examples
    --------
    >>> round(float(zukauskas_nusselt(13943.0, 0.71, 0.0313, 0.0343, "staggered", rows=7, Pr_s=0.70)), 1)
    88.8
    """
    _finite(Re=Re)
    _positive(Pr=Pr, ST=ST, SL=SL)
    if arrangement not in _ARRANGEMENTS:
        raise ValueError(f"arrangement must be one of {list(_ARRANGEMENTS)}")
    _integer("rows", rows, 1)
    Re = np.asarray(Re, dtype=float)
    if np.any((Re < 10) | (Re > 2e6)):
        raise ValueError("Re must lie in [10, 2e6]")
    Pr_ratio = 1.0
    if Pr_s is not None:
        _positive(Pr_s=Pr_s)
        Pr_ratio = np.asarray(Pr, dtype=float) / Pr_s
    upper, C, m = (np.array(column) for column in zip(*_ZUKAUSKAS[arrangement]))
    C = C.copy()
    if arrangement == "staggered":
        C[2] = 0.35 * (ST / SL) ** 0.2 if ST / SL < 2 else 0.40
    band = np.searchsorted(upper, Re)
    C2 = np.interp(min(rows, 20), _ROWS, _ROW_FACTOR[arrangement])
    return C2 * C[band] * Re ** m[band] * np.asarray(Pr, dtype=float) ** 0.36 * Pr_ratio**0.25


_LAMINAR_TUBE = {"flux": 48 / 11, "temperature": 3.6568}


def laminar_tube_nusselt(boundary="flux"):
    r"""Nusselt number of fully developed laminar flow in a circular tube.

    Parameters
    ----------
    boundary : {"flux", "temperature"}, optional
        Uniform wall heat flux (:math:`Nu_D = 48/11 \approx 4.364`) or
        uniform wall temperature (:math:`Nu_D = 3.657`, Graetz's limit).

    Returns
    -------
    float
        Nusselt number :math:`hD/k` [-], based on the mean (bulk) temperature.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Valid for :math:`Re_D \lesssim 2300` beyond
    the thermal entry length (:func:`entry_length`). The uniform-flux value
    follows in closed form from the parabolic profile; the uniform-temperature
    value is twice the square of the first eigenvalue of Graetz's problem
    (scaled by the radius) and has no closed form.

    References
    ----------
    L. Graetz, "Über die Wärmeleitungsfähigkeit von Flüssigkeiten," *Ann.
    Phys. Chem.* 18 (1883) 79-94.

    Examples
    --------
    >>> round(laminar_tube_nusselt("flux"), 3)
    4.364
    """
    if boundary not in _LAMINAR_TUBE:
        raise ValueError(f"boundary must be one of {sorted(_LAMINAR_TUBE)}")
    return _LAMINAR_TUBE[boundary]


def entry_length(Re, D, Pr=1.0):
    r"""Entry length of tube flow: laminar :math:`0.05\,Re_D\,Pr\,D`, turbulent :math:`10D`.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number on the diameter [-], positive; laminar below 2300.
    D : float or array_like
        Tube diameter [m], positive.
    Pr : float or array_like, optional
        Prandtl number [-], positive. The default 1 gives the hydrodynamic
        entry length; the fluid's Prandtl number gives the thermal one.

    Returns
    -------
    float or ndarray
        Distance from the inlet to fully developed conditions [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. The turbulent value is a rule of thumb;
    turbulent profiles develop within 10 to 60 diameters, almost
    independently of the Prandtl number.

    Examples
    --------
    >>> float(entry_length(1000.0, 0.01, 7.0))
    3.5
    """
    _positive(Re=Re, D=D, Pr=Pr)
    Re, D = np.asarray(Re, dtype=float), np.asarray(D, dtype=float)
    return np.where(Re < 2300, 0.05 * Re * np.asarray(Pr, dtype=float) * D, 10 * D)[()]


def hausen_nusselt(Re, Pr, D, L):
    r"""Hausen's average Nusselt number for laminar thermal entry flow in an isothermal tube.

    .. math::

        \bar{Nu}_D = 3.66 + \frac{0.0668\,Gz_D}{1 + 0.04\,Gz_D^{2/3}},
        \qquad Gz_D = \frac{D}{L}\,Re_D\,Pr.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number on the diameter [-], positive; laminar (below 2300).
    Pr : float or array_like
        Prandtl number [-], positive.
    D : float or array_like
        Tube diameter [m], positive.
    L : float or array_like
        Heated length from the inlet [m], positive.

    Returns
    -------
    float or ndarray
        Nusselt number averaged over the length :math:`L` [-]; tends to
        Graetz's 3.66 as :math:`L` grows.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. The velocity profile is fully developed at
    the start of heating (the "thermal entry" or Graetz problem), which also
    fits oils and other high-Prandtl fluids whose flow develops much faster
    than their temperature. Taken from Incropera et al.; Hausen's original
    paper (1943) could not be checked.

    Examples
    --------
    >>> round(float(hausen_nusselt(1000.0, 5.0, 0.01, 1.0)), 3)
    5.825
    """
    _positive(Re=Re, Pr=Pr, D=D, L=L)
    Gz = np.asarray(D, dtype=float) / np.asarray(L, dtype=float) * np.asarray(Re, dtype=float) * np.asarray(Pr, dtype=float)
    return 3.66 + 0.0668 * Gz / (1 + 0.04 * Gz ** (2 / 3))


def lewis_number(alpha, D_AB):
    r"""Lewis number :math:`Le = \alpha/D_{AB} = Sc/Pr`, heat over mass diffusivity.

    Parameters
    ----------
    alpha : float or array_like
        Thermal diffusivity [m^2/s], positive.
    D_AB : float or array_like
        Binary mass diffusivity of species A in B [m^2/s], positive.

    Returns
    -------
    float or ndarray
        Lewis number [-]; about 0.85 for water vapour in air.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(lewis_number(2.2e-5, 2.6e-5)), 3)
    0.846
    """
    _positive(alpha=alpha, D_AB=D_AB)
    return np.asarray(alpha, dtype=float) / D_AB


def chilton_colburn_mass_transfer(h, rho, cp, Le):
    r"""Mass transfer coefficient from a heat transfer coefficient by the Chilton-Colburn analogy.

    Equal Colburn factors :math:`j_H = St\,Pr^{2/3}` and
    :math:`j_m = St_m Sc^{2/3}` give

    .. math::

        \frac{h}{h_m} = \rho\,c_p\,Le^{2/3}.

    Parameters
    ----------
    h : float or array_like
        Heat transfer coefficient [W/(m^2 K)], nonnegative.
    rho : float or array_like
        Density of the fluid [kg/m^3], positive.
    cp : float or array_like
        Specific heat of the fluid [J/(kg K)], positive.
    Le : float or array_like
        Lewis number :math:`\alpha/D_{AB}` [-], positive (see :func:`lewis_number`).

    Returns
    -------
    float or ndarray
        Mass transfer coefficient :math:`h_m` [m/s]; the species flux is
        :math:`h_m(\rho_{A,s} - \rho_{A,\infty})` [kg/(m^2 s)].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Exact between correlations of the form
    :math:`Nu \propto Re^m Pr^{1/3}` and :math:`Sh \propto Re^m Sc^{1/3}`
    (laminar flat plate); a good approximation for turbulent flow with
    :math:`0.6 < Pr, Sc < 60`. Low mass transfer rates are assumed.

    References
    ----------
    T. H. Chilton and A. P. Colburn, "Mass transfer (absorption)
    coefficients: prediction from data on heat transfer and fluid friction,"
    *Ind. Eng. Chem.* 26 (1934) 1183-1187.

    Examples
    --------
    >>> round(float(chilton_colburn_mass_transfer(20.0, 1.16, 1007.0, 0.85)) * 1000, 2)
    19.08
    """
    _nonnegative(h=h)
    _positive(rho=rho, cp=cp, Le=Le)
    return np.asarray(h, dtype=float) / (rho * cp * np.asarray(Le, dtype=float) ** (2 / 3))
