r"""Phase change: moving melt fronts, film condensation and the boiling crisis.

Stefan (1891) analysed the growth of polar sea ice, a conduction problem
whose boundary moves as latent heat is released; Neumann's similarity
solution makes the front advance as :math:`\sqrt t`. Nusselt (1916) solved
laminar film condensation on a vertical wall, and Zuber (1958) explained
the maximum heat flux of Nukiyama's (1934) boiling curve as a
hydrodynamic instability.
"""

from __future__ import annotations

import numpy as np
from scipy import special
from scipy.optimize import brentq

from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.heat_transfer.utils.constants import STANDARD_GRAVITY

__all__ = ["stefan_number", "neumann_stefan_front", "nusselt_film_condensation", "zuber_critical_heat_flux"]


def stefan_number(cp, dT, latent_heat):
    r"""Stefan number :math:`Ste = c_p\,\Delta T/h_{sf}`, sensible over latent heat.

    Parameters
    ----------
    cp : float or array_like
        Specific heat of the phase that conducts [J/(kg K)], positive.
    dT : float or array_like
        Temperature difference driving the phase change [K], nonnegative.
    latent_heat : float or array_like
        Latent heat of fusion (or vaporisation) [J/kg], positive.

    Returns
    -------
    float or ndarray
        Stefan number [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(stefan_number(2100.0, 20.0, 334e3)), 4)
    0.1257
    """
    _positive(cp=cp, latent_heat=latent_heat)
    _nonnegative(dT=dT)
    return cp * np.asarray(dT, dtype=float) / latent_heat


def _neumann_coefficient(Ste):
    """Root of lambda exp(lambda^2) erf(lambda) = Ste / sqrt(pi) (scalar Ste > 0)."""
    target = Ste / np.sqrt(np.pi)
    # The left side increases monotonically from 0; it exceeds Ste/sqrt(pi) at sqrt(Ste) + 1.
    return brentq(lambda lam: lam * np.exp(lam**2) * special.erf(lam) - target, 0.0, np.sqrt(Ste) + 1.0, xtol=1e-15)


def neumann_stefan_front(t, alpha, Ste):
    r"""Position of the phase front in the one-phase Stefan problem (Neumann's solution).

    A liquid at its melting point fills :math:`x > 0`; at :math:`t = 0` the
    face :math:`x = 0` is held :math:`\Delta T` below the melting point and
    a solid layer grows (or, the other way round, a solid at its melting
    point melts). The front is at

    .. math::

        s(t) = 2\lambda\sqrt{\alpha t}, \qquad
        \lambda\,e^{\lambda^2}\operatorname{erf}\lambda = \frac{Ste}{\sqrt\pi},

    and the temperature in the new phase is
    :math:`T_s + \Delta T\,\operatorname{erf}(x/2\sqrt{\alpha t})/\operatorname{erf}\lambda`.

    Parameters
    ----------
    t : float or array_like
        Time [s], nonnegative.
    alpha : float
        Thermal diffusivity of the growing phase [m^2/s], positive.
    Ste : float
        Stefan number of the growing phase [-], positive (see :func:`stefan_number`).

    Returns
    -------
    float or ndarray
        Front position [m]. For small :math:`Ste` this tends to Stefan's
        quasi-steady :math:`\sqrt{2\,Ste\,\alpha t}`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Densities of the two phases are taken
    equal, and the original phase does not conduct.

    References
    ----------
    J. Stefan, "Über die Theorie der Eisbildung, insbesondere über die
    Eisbildung im Polarmeere," *Ann. Phys. Chem.* 42 (1891) 269-286.

    Examples
    --------
    >>> round(float(neumann_stefan_front(86400.0, 1.2e-6, 0.1257)) * 100, 2)
    15.82
    """
    _nonnegative(t=t)
    _positive(alpha=alpha, Ste=Ste)
    lam = _neumann_coefficient(float(Ste))
    return 2 * lam * np.sqrt(alpha * np.asarray(t, dtype=float))


def nusselt_film_condensation(rho_l, rho_v, mu_l, k_l, h_fg, Tsat, Ts, L, g=STANDARD_GRAVITY):
    r"""Nusselt's average coefficient for laminar film condensation on a vertical wall.

    .. math::

        \bar h = 0.943\left[\frac{g\,\rho_l(\rho_l - \rho_v)\,k_l^3\,h_{fg}}
                 {\mu_l\,L\,(T_{sat} - T_s)}\right]^{1/4}

    The film thickness grows as
    :math:`\delta(x) = [4\mu_l k_l (T_{sat}-T_s)\,x/(g\rho_l(\rho_l-\rho_v)h_{fg})]^{1/4}`,
    and the local coefficient is :math:`k_l/\delta`.

    Parameters
    ----------
    rho_l, rho_v : float
        Liquid and vapour densities [kg/m^3], with ``rho_l > rho_v >= 0``.
    mu_l : float
        Liquid dynamic viscosity [Pa s], positive.
    k_l : float
        Liquid conductivity [W/(m K)], positive.
    h_fg : float
        Latent heat [J/kg], positive. Rohsenow's
        :math:`h_{fg} + 0.68\,c_{p,l}(T_{sat}-T_s)` accounts for subcooling of the film.
    Tsat : float
        Saturation temperature [K].
    Ts : float or array_like
        Wall temperature [K], below ``Tsat``.
    L : float or array_like
        Wall height [m], positive.
    g : float, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Average heat-transfer coefficient over the height [W/(m^2 K)].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Laminar film (film Reynolds number
    below about 1800); ripples above about 30 raise the coefficient by
    10-20%. Properties at the film temperature.

    References
    ----------
    W. Nusselt, "Die Oberflächenkondensation des Wasserdampfes," *Z. VDI*
    60 (1916) 541-546, 569-575.

    Examples
    --------
    >>> round(float(nusselt_film_condensation(960.0, 0.6, 2.8e-4, 0.68, 2.257e6, 373.15, 363.15, 0.5)))
    7757
    """
    _positive(rho_l=rho_l, mu_l=mu_l, k_l=k_l, h_fg=h_fg, Tsat=Tsat, L=L, g=g)
    _nonnegative(rho_v=rho_v, Ts=Ts)
    if np.any(np.asarray(rho_v) >= rho_l):
        raise ValueError("rho_l must exceed rho_v")
    dT = Tsat - np.asarray(Ts, dtype=float)
    if np.any(dT <= 0):
        raise ValueError("Ts must be below Tsat")
    return 0.943 * (g * rho_l * (rho_l - rho_v) * k_l**3 * h_fg / (mu_l * np.asarray(L, dtype=float) * dT)) ** 0.25


def zuber_critical_heat_flux(rho_l, rho_v, h_fg, sigma, g=STANDARD_GRAVITY, C=np.pi / 24):
    r"""Zuber's critical (maximum) heat flux in saturated pool boiling.

    .. math::

        q''_{max} = C\,h_{fg}\,\rho_v\left[\frac{\sigma g(\rho_l - \rho_v)}{\rho_v^2}\right]^{1/4}

    Parameters
    ----------
    rho_l, rho_v : float or array_like
        Liquid and vapour densities [kg/m^3], with ``rho_l > rho_v > 0``.
    h_fg : float or array_like
        Latent heat of vaporisation [J/kg], positive.
    sigma : float or array_like
        Surface tension [N/m], positive.
    g : float, optional
        Gravitational acceleration [m/s^2], positive.
    C : float, optional
        Constant [-]: Zuber's :math:`\pi/24 \approx 0.131`; Lienhard and Dhir's
        0.149 fits large flat heaters better.

    Returns
    -------
    float or ndarray
        Critical heat flux [W/m^2]; beyond it a vapour film blankets the
        heater and its temperature jumps (burnout).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Large horizontal heater, saturated liquid.

    References
    ----------
    S. Nukiyama, "The maximum and minimum values of the heat Q transmitted
    from metal to boiling water under atmospheric pressure," *J. Soc. Mech.
    Eng. Japan* 37 (1934) 367-374 (transl. *Int. J. Heat Mass Transfer* 9
    (1966) 1419-1433); N. Zuber, "On the stability of boiling heat
    transfer," *Trans. ASME* 80 (1958) 711-720.

    Examples
    --------
    >>> round(float(zuber_critical_heat_flux(957.9, 0.5955, 2.257e6, 0.0589)) / 1e6, 3)
    1.106
    """
    _positive(rho_l=rho_l, rho_v=rho_v, h_fg=h_fg, sigma=sigma, g=g, C=C)
    rho_l, rho_v = np.asarray(rho_l, dtype=float), np.asarray(rho_v, dtype=float)
    if np.any(rho_v >= rho_l):
        raise ValueError("rho_l must exceed rho_v")
    return C * h_fg * rho_v * (sigma * g * (rho_l - rho_v) / rho_v**2) ** 0.25
