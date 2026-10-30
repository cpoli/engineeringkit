r"""Thermal radiation: blackbody laws, gray-surface exchange and view factors.

Kirchhoff (1859) proved that a body's emissivity equals its absorptivity,
which makes the blackbody the universal reference. Stefan (1879) and
Boltzmann (1884) found its total emissive power :math:`\sigma T^4`; Wien
(1893) showed its spectrum shifts as :math:`\lambda_{max}T = b`; Planck
(1900) gave its exact shape. Nusselt (1928) reduced the geometry of
exchange between surfaces to a projection on a unit sphere, and Oppenheim
(1956) the exchange between gray surfaces to an electrical network of
surface and space resistances.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.heat_transfer.core.base import EnclosureResult
from engineeringkit.heat_transfer.utils.constants import BOLTZMANN, PLANCK, SPEED_OF_LIGHT, STEFAN_BOLTZMANN, WIEN_DISPLACEMENT

__all__ = [
    "radiation_heat_rate",
    "radiation_coefficient",
    "stefan_boltzmann_emissive_power",
    "blackbody_peak_wavelength",
    "planck_spectral_emissive_power",
    "view_factor_parallel_plates",
    "view_factor_coaxial_disks",
    "view_factor_reciprocity",
    "two_surface_radiation",
    "radiation_shield_heat_flux",
    "radiation_enclosure",
]


def radiation_heat_rate(epsilon, A, Ts, Tsur, sigma=STEFAN_BOLTZMANN):
    r"""Net radiation from a small gray surface to large surroundings, :math:`\dot Q = \varepsilon\sigma A(T_s^4 - T_{sur}^4)`.

    Parameters
    ----------
    epsilon : float or array_like
        Surface emissivity [-], in [0, 1]; by Kirchhoff's law also its absorptivity.
    A : float or array_like
        Surface area [m^2], positive.
    Ts : float or array_like
        Absolute surface temperature [K], nonnegative.
    Tsur : float or array_like
        Absolute temperature of the surroundings [K], nonnegative.
    sigma : float, optional
        Stefan-Boltzmann constant [W/(m^2 K^4)].

    Returns
    -------
    float or ndarray
        Net heat rate leaving the surface [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. The surroundings are a large isothermal
    enclosure; Celsius temperatures give wrong answers silently.

    References
    ----------
    G. Kirchhoff, "Über das Verhältniss zwischen dem Emissionsvermögen und
    dem Absorptionsvermögen der Körper für Wärme und Licht," *Ann. Phys.
    Chem.* 109 (1860) 275-301.

    Examples
    --------
    >>> round(float(radiation_heat_rate(0.8, 1.0, 400.0, 300.0)), 1)
    793.9
    """
    _finite(epsilon=epsilon, A=A, Ts=Ts, Tsur=Tsur, sigma=sigma)
    _interval("epsilon", epsilon, 0, 1)
    _positive(A=A, sigma=sigma)
    _nonnegative(Ts=Ts, Tsur=Tsur)
    return epsilon * sigma * A * (np.asarray(Ts, dtype=float) ** 4 - np.asarray(Tsur, dtype=float) ** 4)


def radiation_coefficient(epsilon, Ts, Tsur, sigma=STEFAN_BOLTZMANN):
    r"""Radiation heat-transfer coefficient :math:`h_r = \varepsilon\sigma(T_s + T_{sur})(T_s^2 + T_{sur}^2)`.

    With it, radiation can be written like convection,
    :math:`\dot Q = h_r A (T_s - T_{sur})`, and added to :math:`h` in parallel.

    Parameters
    ----------
    epsilon : float or array_like
        Surface emissivity [-], in [0, 1].
    Ts : float or array_like
        Absolute surface temperature [K], nonnegative.
    Tsur : float or array_like
        Absolute temperature of the surroundings [K], nonnegative.
    sigma : float, optional
        Stefan-Boltzmann constant [W/(m^2 K^4)].

    Returns
    -------
    float or ndarray
        Radiation coefficient [W/(m^2 K)].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(radiation_coefficient(0.8, 400.0, 300.0)), 3)
    7.939
    """
    _finite(epsilon=epsilon, Ts=Ts, Tsur=Tsur, sigma=sigma)
    _interval("epsilon", epsilon, 0, 1)
    _positive(sigma=sigma)
    _nonnegative(Ts=Ts, Tsur=Tsur)
    Ts, Tsur = np.asarray(Ts, dtype=float), np.asarray(Tsur, dtype=float)
    return epsilon * sigma * (Ts + Tsur) * (Ts**2 + Tsur**2)


def stefan_boltzmann_emissive_power(T, sigma=STEFAN_BOLTZMANN):
    r"""Total emissive power of a blackbody, :math:`E_b = \sigma T^4` (Stefan-Boltzmann law).

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], nonnegative.
    sigma : float, optional
        Stefan-Boltzmann constant [W/(m^2 K^4)].

    Returns
    -------
    float or ndarray
        Emissive power [W/m^2].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    J. Stefan, "Über die Beziehung zwischen der Wärmestrahlung und der
    Temperatur," *Sitzungsber. Akad. Wiss. Wien* 79 (1879) 391-428;
    L. Boltzmann, "Ableitung des Stefan'schen Gesetzes," *Ann. Phys. Chem.*
    22 (1884) 291-294.

    Examples
    --------
    >>> round(float(stefan_boltzmann_emissive_power(5772.0)) / 1e6, 2)
    62.94
    """
    _finite(T=T, sigma=sigma)
    _nonnegative(T=T)
    _positive(sigma=sigma)
    return sigma * np.asarray(T, dtype=float) ** 4


def blackbody_peak_wavelength(T):
    r"""Wavelength of maximum blackbody spectral emissive power, :math:`\lambda_{max} = b/T` (Wien's displacement law).

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], positive.

    Returns
    -------
    float or ndarray
        Peak wavelength [m]; :math:`b = 2.898\times10^{-3}` m K.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    W. Wien, "Eine neue Beziehung der Strahlung schwarzer Körper zum zweiten
    Hauptsatz der Wärmetheorie," *Sitzungsber. Preuss. Akad. Wiss. Berlin*
    (1893) 55-62.

    Examples
    --------
    >>> round(float(blackbody_peak_wavelength(5772.0)) * 1e9, 1)
    502.0
    """
    _finite(T=T)
    _positive(T=T)
    return WIEN_DISPLACEMENT / np.asarray(T, dtype=float)


def planck_spectral_emissive_power(wavelength, T):
    r"""Planck's spectral emissive power of a blackbody.

    .. math::

        E_{\lambda,b} = \frac{2\pi h c^2}{\lambda^5\,[\exp(hc/\lambda k T) - 1]}

    Parameters
    ----------
    wavelength : float or array_like
        Wavelength [m], positive.
    T : float or array_like
        Absolute temperature [K], positive.

    Returns
    -------
    float or ndarray
        Spectral emissive power [W/(m^2 m)], broadcast over the inputs. Its
        integral over wavelength is :math:`\sigma T^4`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Written with :math:`e^{-x}/(1 - e^{-x})`
    so that very short wavelengths underflow to zero instead of overflowing.

    References
    ----------
    M. Planck, "Über das Gesetz der Energieverteilung im Normalspectrum,"
    *Ann. Phys.* 4 (1901) 553-563.

    Examples
    --------
    >>> f"{float(planck_spectral_emissive_power(10e-6, 300.0)):.4e}"
    '3.1177e+07'
    """
    _positive(wavelength=wavelength, T=T)
    lam, T = np.asarray(wavelength, dtype=float), np.asarray(T, dtype=float)
    x = PLANCK * SPEED_OF_LIGHT / (lam * BOLTZMANN * T)
    return 2 * np.pi * PLANCK * SPEED_OF_LIGHT**2 / lam**5 * np.exp(-x) / -np.expm1(-x)


def view_factor_parallel_plates(a, b, L):
    r"""View factor between two aligned, directly opposed parallel rectangles.

    With :math:`X = a/L` and :math:`Y = b/L`,

    .. math::

        F_{12} = \frac{2}{\pi XY}\Big\{\ln\Big[\frac{(1+X^2)(1+Y^2)}{1+X^2+Y^2}\Big]^{1/2}
        + X\sqrt{1+Y^2}\arctan\frac{X}{\sqrt{1+Y^2}}
        + Y\sqrt{1+X^2}\arctan\frac{Y}{\sqrt{1+X^2}} - X\arctan X - Y\arctan Y\Big\}.

    Parameters
    ----------
    a, b : float or array_like
        Side lengths of each rectangle [m], positive.
    L : float or array_like
        Separation [m], positive.

    Returns
    -------
    float or ndarray
        Fraction of the diffuse radiation leaving one rectangle that strikes the other [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    W. Nusselt, "Graphische Bestimmung des Winkelverhältnisses bei der
    Wärmestrahlung," *Z. VDI* 72 (1928) 673; J. R. Howell, *A Catalog of
    Radiation Heat Transfer Configuration Factors*, case C-11.

    Examples
    --------
    >>> round(float(view_factor_parallel_plates(1.0, 1.0, 1.0)), 4)
    0.1998
    """
    _positive(a=a, b=b, L=L)
    X, Y = np.asarray(a, dtype=float) / L, np.asarray(b, dtype=float) / L
    sx, sy = np.sqrt(1 + X**2), np.sqrt(1 + Y**2)
    bracket = (
        0.5 * np.log(sx**2 * sy**2 / (1 + X**2 + Y**2))
        + X * sy * np.arctan(X / sy)
        + Y * sx * np.arctan(Y / sx)
        - X * np.arctan(X)
        - Y * np.arctan(Y)
    )
    return 2 / (np.pi * X * Y) * bracket


def view_factor_coaxial_disks(r1, r2, L):
    r"""View factor from one disk to a parallel, coaxial disk.

    With :math:`R_1 = r_1/L`, :math:`R_2 = r_2/L` and
    :math:`S = 1 + (1 + R_2^2)/R_1^2`,
    :math:`F_{12} = \tfrac12\{S - [S^2 - 4(r_2/r_1)^2]^{1/2}\}`.

    Parameters
    ----------
    r1 : float or array_like
        Radius of the emitting disk [m], positive.
    r2 : float or array_like
        Radius of the receiving disk [m], positive.
    L : float or array_like
        Separation [m], positive.

    Returns
    -------
    float or ndarray
        View factor :math:`F_{12}` [-]; reciprocity gives :math:`r_1^2F_{12} = r_2^2F_{21}`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    W. Nusselt, "Graphische Bestimmung des Winkelverhältnisses bei der
    Wärmestrahlung," *Z. VDI* 72 (1928) 673; J. R. Howell, *A Catalog of
    Radiation Heat Transfer Configuration Factors*, case C-41.

    Examples
    --------
    >>> round(float(view_factor_coaxial_disks(1.0, 1.0, 1.0)), 4)
    0.382
    """
    _positive(r1=r1, r2=r2, L=L)
    r1, r2 = np.asarray(r1, dtype=float), np.asarray(r2, dtype=float)
    R1, R2 = r1 / L, r2 / L
    S = 1 + (1 + R2**2) / R1**2
    return 0.5 * (S - np.sqrt(S**2 - 4 * (r2 / r1) ** 2))


def view_factor_reciprocity(F_ij, A_i, A_j):
    r"""View factor :math:`F_{ji}` from its reverse by reciprocity, :math:`A_iF_{ij} = A_jF_{ji}`.

    Parameters
    ----------
    F_ij : float or array_like
        View factor from surface :math:`i` to surface :math:`j` [-], in [0, 1].
    A_i, A_j : float or array_like
        Areas of the two surfaces [m^2], positive.

    Returns
    -------
    float or ndarray
        View factor :math:`F_{ji}` [-].

    Raises
    ------
    ValueError
        If the result exceeds one, which no pair of real surfaces allows.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. With the summation rule
    :math:`\sum_j F_{ij} = 1` over an enclosure, reciprocity gives every view
    factor of a two- or three-surface enclosure from its areas.

    Examples
    --------
    >>> float(view_factor_reciprocity(1.0, 1.0, 4.0))
    0.25
    """
    _interval("F_ij", F_ij, 0, 1)
    _positive(A_i=A_i, A_j=A_j)
    F_ji = A_i * np.asarray(F_ij, dtype=float) / A_j
    if np.any(F_ji > 1 + 1e-12):
        raise ValueError("reciprocity gives F_ji > 1: check F_ij and the areas")
    return F_ji[()]


def two_surface_radiation(T1, T2, eps1, eps2, A1, A2, F12=1.0, sigma=STEFAN_BOLTZMANN):
    r"""Net radiation exchange in a two-surface gray enclosure.

    .. math::

        \dot Q_{12} = \frac{\sigma(T_1^4 - T_2^4)}
        {\dfrac{1 - \varepsilon_1}{\varepsilon_1A_1} + \dfrac{1}{A_1F_{12}} + \dfrac{1 - \varepsilon_2}{\varepsilon_2A_2}}

    Parameters
    ----------
    T1, T2 : float or array_like
        Absolute temperatures of the two surfaces [K], nonnegative.
    eps1, eps2 : float
        Emissivities [-], in (0, 1].
    A1, A2 : float
        Areas [m^2], positive.
    F12 : float, optional
        View factor from surface 1 to surface 2 [-], in (0, 1].
    sigma : float, optional
        Stefan-Boltzmann constant [W/(m^2 K^4)].

    Returns
    -------
    float or ndarray
        Net heat rate from surface 1 to surface 2 [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Diffuse gray surfaces, each isothermal
    with uniform radiosity. Large parallel plates have :math:`A_1 = A_2` and
    :math:`F_{12} = 1`; a convex body inside a large enclosure
    (:math:`A_1/A_2 \to 0`) recovers :func:`radiation_heat_rate`.

    References
    ----------
    A. K. Oppenheim, "Radiation analysis by the network method," *Trans.
    ASME* 78 (1956) 725-735.

    Examples
    --------
    >>> round(float(two_surface_radiation(500.0, 300.0, 0.8, 0.8, 1.0, 1.0)), 1)
    2056.5
    """
    _finite(T1=T1, T2=T2, sigma=sigma)
    _nonnegative(T1=T1, T2=T2)
    _interval("eps1", eps1, 0, 1, closed="right")
    _interval("eps2", eps2, 0, 1, closed="right")
    _interval("F12", F12, 0, 1, closed="right")
    _positive(A1=A1, A2=A2, sigma=sigma)
    resistance = (1 - eps1) / (eps1 * A1) + 1 / (A1 * F12) + (1 - eps2) / (eps2 * A2)
    return sigma * (np.asarray(T1, dtype=float) ** 4 - np.asarray(T2, dtype=float) ** 4) / resistance


def radiation_shield_heat_flux(T1, T2, eps1, eps2, shields=(), sigma=STEFAN_BOLTZMANN):
    r"""Radiation flux between large parallel plates separated by thin radiation shields.

    Each gap adds a resistance :math:`1/\varepsilon_a + 1/\varepsilon_b - 1`
    per unit area, so with shields whose faces have emissivities
    :math:`\varepsilon_{k,a}` and :math:`\varepsilon_{k,b}`

    .. math::

        q'' = \frac{\sigma(T_1^4 - T_2^4)}{\dfrac{1}{\varepsilon_1} + \dfrac{1}{\varepsilon_2} - 1
              + \sum_k\left(\dfrac{1}{\varepsilon_{k,a}} + \dfrac{1}{\varepsilon_{k,b}} - 1\right)}.

    Parameters
    ----------
    T1, T2 : float or array_like
        Absolute temperatures of the two plates [K], nonnegative.
    eps1, eps2 : float
        Plate emissivities [-], in (0, 1].
    shields : sequence of float or (float, float), optional
        One entry per shield: a single emissivity for both faces, or the pair
        ``(eps_a, eps_b)`` of the faces towards plate 1 and plate 2, each in
        (0, 1]. Empty means no shield.
    sigma : float, optional
        Stefan-Boltzmann constant [W/(m^2 K^4)].

    Returns
    -------
    float or ndarray
        Net heat flux from plate 1 to plate 2 [W/m^2]. :math:`N` shields with
        the plates' emissivity divide it by :math:`N + 1`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Shields are thin, conduct well across their
    thickness and fill the gap (all view factors one). Shield temperatures
    follow from the flux, one resistance at a time.

    References
    ----------
    A. K. Oppenheim, "Radiation analysis by the network method," *Trans.
    ASME* 78 (1956) 725-735.

    Examples
    --------
    >>> round(float(radiation_shield_heat_flux(500.0, 300.0, 0.8, 0.8, [0.05])), 1)
    76.2
    """
    _finite(T1=T1, T2=T2, sigma=sigma)
    _nonnegative(T1=T1, T2=T2)
    _positive(sigma=sigma)
    faces = [eps1, eps2]
    for shield in shields:
        faces += list(shield) if np.ndim(shield) else [shield, shield]
    for eps in faces:
        _interval("emissivity", eps, 0, 1, closed="right")
    resistance = sum(1 / eps for eps in faces) - (1 + len(shields))
    return sigma * (np.asarray(T1, dtype=float) ** 4 - np.asarray(T2, dtype=float) ** 4) / resistance


def radiation_enclosure(T, epsilon, A, F, reradiating=None, sigma=STEFAN_BOLTZMANN):
    r"""Net radiation exchange in an enclosure of diffuse gray surfaces (Oppenheim's network).

    Each surface's radiosity :math:`J_i` is linked to its blackbody emissive
    power by the surface resistance :math:`(1 - \varepsilon_i)/(\varepsilon_iA_i)`
    and to every other radiosity by the space resistance :math:`1/(A_iF_{ij})`.
    The node equations reduce to

    .. math::

        J_i - (1 - \varepsilon_i)\sum_j F_{ij}J_j = \varepsilon_i\sigma T_i^4,

    with :math:`\varepsilon_i` replaced by 0 for a reradiating (insulated)
    surface, and the net radiation leaving surface :math:`i` is
    :math:`\dot Q_i = A_i(J_i - \sum_j F_{ij}J_j)`.

    Parameters
    ----------
    T : array_like
        Absolute temperature of each surface [K], nonnegative; ignored for
        reradiating surfaces. ``n`` surfaces, at least two.
    epsilon : array_like
        Emissivity of each surface [-], in [0, 1]; 1 is black.
    A : array_like
        Area of each surface [m^2], positive.
    F : array_like
        View-factor matrix [-], shape ``(n, n)``; ``F[i, j]`` is :math:`F_{ij}`.
        Each row must sum to one (the enclosure is closed) and
        :math:`A_iF_{ij} = A_jF_{ji}` must hold, both to 1e-6.
    reradiating : array_like of bool, optional
        Surfaces that are insulated (net radiation zero); their temperature
        is solved for.
    sigma : float, optional
        Stefan-Boltzmann constant [W/(m^2 K^4)].

    Returns
    -------
    EnclosureResult
        Radiosities [W/m^2], net heat rates [W] summing to zero, and the
        temperatures [K] with the reradiating ones filled in.

    Raises
    ------
    ValueError
        If the shapes disagree, a view factor lies outside [0, 1], or the
        summation or reciprocity rule fails.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Each surface is isothermal, diffuse and
    gray with uniform radiosity; the medium between them does not
    participate.

    References
    ----------
    A. K. Oppenheim, "Radiation analysis by the network method," *Trans.
    ASME* 78 (1956) 725-735.

    Examples
    --------
    >>> F = [[0.0, 1.0], [0.25, 0.75]]
    >>> result = radiation_enclosure([500.0, 300.0], [0.8, 0.8], [1.0, 4.0], F)
    >>> round(float(result.heat_rate[0]), 1)
    2350.2
    """
    T = _vector("T", T, min_length=2)
    n = T.size
    epsilon, A = _vector("epsilon", epsilon), _vector("A", A)
    F = np.asarray(F, dtype=float)
    if epsilon.size != n or A.size != n or F.shape != (n, n):
        raise ValueError("T, epsilon and A need one entry per surface and F shape (n, n)")
    _nonnegative(T=T)
    _interval("epsilon", epsilon, 0, 1)
    _positive(A=A, sigma=sigma)
    _interval("F", F, 0, 1)
    if not np.allclose(F.sum(axis=1), 1.0, rtol=0, atol=1e-6):
        raise ValueError("each row of F must sum to one (summation rule)")
    exchange = A[:, None] * F
    if not np.allclose(exchange, exchange.T, rtol=0, atol=1e-6 * A.max()):
        raise ValueError("F must satisfy A_i F_ij = A_j F_ji (reciprocity)")
    insulated = np.zeros(n, dtype=bool) if reradiating is None else np.asarray(reradiating, dtype=bool)
    if insulated.shape != (n,):
        raise ValueError("reradiating needs one entry per surface")
    eps = np.where(insulated, 0.0, epsilon)
    J = np.linalg.solve(np.eye(n) - (1 - eps)[:, None] * F, eps * sigma * T**4)
    heat_rate = A * (J - F @ J)
    temperature = np.where(insulated, (np.maximum(J, 0.0) / sigma) ** 0.25, T)
    return EnclosureResult(J, heat_rate, temperature)
