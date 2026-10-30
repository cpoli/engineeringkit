r"""Maxwell's relations and the general property relations that follow from them.

Because :math:`u`, :math:`h`, the Helmholtz function :math:`a` and the Gibbs
function :math:`g` are properties, their mixed second derivatives are equal.
Maxwell wrote the four consequences,

.. math:: \left(\frac{\partial T}{\partial v}\right)_s = -\left(\frac{\partial p}{\partial s}\right)_v,\quad
          \left(\frac{\partial T}{\partial p}\right)_s = \left(\frac{\partial v}{\partial s}\right)_p,\quad
          \left(\frac{\partial s}{\partial v}\right)_T = \left(\frac{\partial p}{\partial T}\right)_v,\quad
          \left(\frac{\partial s}{\partial p}\right)_T = -\left(\frac{\partial v}{\partial T}\right)_p,

which replace unmeasurable entropy derivatives by measurable :math:`p`-:math:`v`-:math:`T`
ones. Two of their best-known results are the difference of the specific heats
and the general Joule-Thomson coefficient. Both need only the volume expansivity
:math:`\beta = (1/v)(\partial v/\partial T)_p` and the isothermal
compressibility :math:`\kappa_T = -(1/v)(\partial v/\partial p)_T`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive

__all__ = ["cp_minus_cv", "joule_thomson_coefficient"]


def cp_minus_cv(T, v, beta, kappa_T):
    r"""Difference of the specific heats, :math:`c_p - c_v = Tv\beta^2/\kappa_T`.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], positive.
    v : float or array_like
        Specific volume [m^3/kg], positive.
    beta : float or array_like
        Volume expansivity :math:`(1/v)(\partial v/\partial T)_p` [1/K].
    kappa_T : float or array_like
        Isothermal compressibility :math:`-(1/v)(\partial v/\partial p)_T`
        [1/Pa], positive.

    Returns
    -------
    float or ndarray
        :math:`c_p - c_v` [J/(kg K)]: never negative, the gas constant
        :math:`R` for an ideal gas, and zero for water at 4 °C where
        :math:`\beta = 0`.

    Notes
    -----
    Follows from :math:`(\partial s/\partial v)_T = (\partial p/\partial T)_v`
    and :math:`(\partial s/\partial p)_T = -(\partial v/\partial T)_p`.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    J. C. Maxwell, *Theory of Heat*, London: Longmans, Green (1871), Ch. 9.

    Examples
    --------
    >>> round(float(cp_minus_cv(300.0, 287.0 * 300.0 / 1e5, 1 / 300.0, 1 / 1e5)), 6)
    287.0
    """
    _positive(T=T, v=v, kappa_T=kappa_T)
    _finite(beta=beta)
    return (np.asarray(T, dtype=float) * v * np.asarray(beta, dtype=float) ** 2 / kappa_T)[()]


def joule_thomson_coefficient(T, v, beta, cp):
    r"""Joule-Thomson coefficient of any fluid, :math:`\mu_{JT} = v(T\beta - 1)/c_p`.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], positive.
    v : float or array_like
        Specific volume [m^3/kg], positive.
    beta : float or array_like
        Volume expansivity [1/K].
    cp : float or array_like
        Isobaric specific heat [J/(kg K)], positive.

    Returns
    -------
    float or ndarray
        :math:`(\partial T/\partial p)_h` [K/Pa]: zero for an ideal gas
        (:math:`T\beta = 1`), positive where throttling cools.

    Notes
    -----
    From :math:`dh = c_p\,dT + [v - T(\partial v/\partial T)_p]\,dp`, whose
    second term comes from the Maxwell relation
    :math:`(\partial s/\partial p)_T = -(\partial v/\partial T)_p`. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    J. C. Maxwell (1871), as above.

    Examples
    --------
    >>> float(joule_thomson_coefficient(300.0, 0.861, 1 / 300.0, 1005.0))
    0.0
    """
    _positive(T=T, v=v, cp=cp)
    _finite(beta=beta)
    return (np.asarray(v, dtype=float) * (np.asarray(T, dtype=float) * beta - 1) / cp)[()]
