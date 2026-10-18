r"""Vapour pressure: the Clausius-Clapeyron relation and Antoine's equation.

Clapeyron related the slope of the boiling curve to the latent heat; Clausius
simplified it for an ideal vapour and a negligible liquid volume, so that
:math:`\ln p` is linear in :math:`1/T`. Antoine's empirical third constant
bends that line to fit measured vapour pressures over a wide range.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive

__all__ = ["clausius_clapeyron_pressure", "antoine_vapor_pressure"]


def clausius_clapeyron_pressure(T, T_ref, p_ref, h_fg, R_gas):
    r"""Vapour pressure from the integrated Clausius-Clapeyron relation.

    .. math:: p = p_{ref}\exp\left[-\frac{h_{fg}}{R}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right]

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], positive.
    T_ref : float
        Absolute temperature of a known saturation state [K], positive.
    p_ref : float
        Saturation pressure at ``T_ref`` [Pa], positive.
    h_fg : float
        Latent heat of vaporization [J/kg], positive, taken as constant.
    R_gas : float
        Specific gas constant of the vapour [J/(kg K)], positive (461.5 for
        water).

    Returns
    -------
    float or ndarray
        Saturation pressure [Pa].

    Notes
    -----
    Assumes the vapour is an ideal gas, the liquid volume is negligible and
    :math:`h_{fg}` does not vary with temperature. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    É. Clapeyron, "Mémoire sur la puissance motrice de la chaleur," *J. École
    Polytechnique* 14 (1834) 153-190; R. Clausius, "Über die bewegende Kraft
    der Wärme," *Ann. Phys.* 79 (1850) 368-397, 500-524.

    Examples
    --------
    >>> round(float(clausius_clapeyron_pressure(363.15, 373.15, 101325.0, 2.257e6, 461.5)), -2)
    70600.0
    """
    _positive(T=T, T_ref=T_ref, p_ref=p_ref, h_fg=h_fg, R_gas=R_gas)
    return p_ref * np.exp(-h_fg / R_gas * (1 / np.asarray(T, dtype=float) - 1 / T_ref))


def antoine_vapor_pressure(T, A, B, C):
    r"""Antoine's vapour-pressure equation, :math:`\log_{10} p = A - B/(C + T)`.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], with :math:`C + T > 0`.
    A : float
        Constant [-] for ``p`` in Pa.
    B : float
        Constant [K].
    C : float
        Constant [K] for ``T`` in kelvin.

    Returns
    -------
    float or ndarray
        Vapour pressure [Pa].

    Notes
    -----
    Tabulated constants are often for mmHg and °C. To convert, add
    :math:`\log_{10} 133.322 = 2.12490` to :math:`A` and subtract 273.15
    from :math:`C`; :math:`B` is unchanged. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    C. Antoine, "Tensions des vapeurs; nouvelle relation entre les tensions
    et les températures," *C. R. Acad. Sci.* 107 (1888) 681-684, 778-780,
    836-837.

    Examples
    --------
    >>> round(float(antoine_vapor_pressure(373.15, 8.07131 + 2.12490, 1730.63, 233.426 - 273.15)), -1)
    101340.0
    """
    _finite(A=A, B=B, C=C)
    _positive(T=T)
    shifted = C + np.asarray(T, dtype=float)
    _positive(**{"C + T": shifted})
    return 10.0 ** (A - B / shifted)
