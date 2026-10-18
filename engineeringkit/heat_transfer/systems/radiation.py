r"""Thermal radiation: blackbody laws, gray-surface exchange and view factors.

Kirchhoff (1859) proved that a body's emissivity equals its absorptivity,
which makes the blackbody the universal reference. Stefan (1879) and
Boltzmann (1884) found its total emissive power :math:`\sigma T^4`; Wien
(1893) showed its spectrum shifts as :math:`\lambda_{max}T = b`; Planck
(1900) gave its exact shape. Nusselt (1928) reduced the geometry of
exchange between surfaces to a projection on a unit sphere.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.heat_transfer.utils.constants import BOLTZMANN, PLANCK, SPEED_OF_LIGHT, STEFAN_BOLTZMANN, WIEN_DISPLACEMENT

__all__ = [
    "radiation_heat_rate",
    "radiation_coefficient",
    "stefan_boltzmann_emissive_power",
    "blackbody_peak_wavelength",
    "planck_spectral_emissive_power",
    "view_factor_parallel_plates",
    "view_factor_coaxial_disks",
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
