r"""The Peng-Robinson cubic equation of state and refrigerant properties from it.

Peng and Robinson kept van der Waals's form, a repulsive term and an
attractive term, but chose the attraction so that one equation fits both the
vapour pressure and the liquid density of hydrocarbons:

.. math:: p = \frac{RT}{v - b} - \frac{a\,\alpha(T)}{v^2 + 2bv - b^2},
   \qquad \alpha = [1 + \kappa(1 - \sqrt{T/T_c})]^2,

with :math:`a` and :math:`b` fixed by the critical point and :math:`\kappa`
by Pitzer's acentric factor. Saturation follows from equal fugacities of
liquid and vapour, and enthalpy and entropy from the ideal gas plus departure
functions. Three constants and an ideal-gas heat capacity are enough to run
a refrigeration cycle on any fluid, with no property tables.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit.thermodynamics.core.base import CubicFluid, RefrigerantState
from engineeringkit.thermodynamics.systems.gases import R_UNIVERSAL

__all__ = [
    "R134A",
    "peng_robinson_pressure",
    "peng_robinson_compressibility",
    "peng_robinson_saturation_pressure",
    "refrigerant_properties",
]

R134A = CubicFluid(name="R-134a", Tc=374.2, pc=4.059e6, omega=0.327, molar_mass=0.10203, cp0=852.0)
"""R-134a (1,1,1,2-tetrafluoroethane). Critical point, molar mass and
ideal-gas :math:`c_p` at 300 K from Çengel and Boles, Tables A-1 and A-2a
(the NIST WebBook averages are 374.1 ± 0.2 K and 4.06 ± 0.02 MPa); acentric
factor from Poling, Prausnitz and O'Connell, *The Properties of Gases and
Liquids*, 5th ed., Appendix A."""

_SQRT2 = math.sqrt(2.0)
_T_LOW = 0.45  # lowest reduced temperature handled (R-134a's triple point is 0.454 Tc)
_T_HIGH = 2.5  # highest reduced temperature searched when inverting s or h
_T_TOP = 0.9999  # highest reduced saturation temperature: the spinodals merge at the critical point


def _constants(T, Tc, pc, omega):
    """Molar :math:`a\\alpha(T)` [Pa m^6/mol^2], :math:`b` [m^3/mol] and :math:`d(a\\alpha)/dT`."""
    a = 0.45724 * (R_UNIVERSAL * Tc) ** 2 / pc
    b = 0.07780 * R_UNIVERSAL * Tc / pc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    root = 1 + kappa * (1 - np.sqrt(T / Tc))
    return a * root**2, b, -a * kappa * root / np.sqrt(T * Tc)


def _check(T, Tc, pc, omega):
    _positive(T=T, Tc=Tc, pc=pc)
    _interval("omega", omega, -0.5, 1.5)


def peng_robinson_pressure(T, v, Tc, pc, omega):
    r"""Peng-Robinson pressure from temperature and molar volume.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], positive.
    v : float or array_like
        Molar volume [m^3/mol], larger than the co-volume
        :math:`b = 0.07780\,RT_c/p_c`.
    Tc : float
        Critical temperature [K], positive.
    pc : float
        Critical pressure [Pa], positive.
    omega : float
        Acentric factor [-], in [-0.5, 1.5].

    Returns
    -------
    float or ndarray
        Pressure [Pa].

    Notes
    -----
    At the critical point the critical isotherm has a horizontal
    inflection and the critical compressibility is 0.3074, nearer real
    fluids than van der Waals's 3/8. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    D.-Y. Peng and D. B. Robinson, "A new two-constant equation of state,"
    *Ind. Eng. Chem. Fundam.* 15 (1976) 59-64.

    Examples
    --------
    >>> round(float(peng_robinson_pressure(300.0, 0.024, 190.56, 4.599e6, 0.011)), -1)
    103700.0
    """
    _check(T, Tc, pc, omega)
    _positive(v=v)
    a_alpha, b, _ = _constants(np.asarray(T, dtype=float), Tc, pc, omega)
    v = np.asarray(v, dtype=float)
    _positive(**{"v - b": v - b})
    return (R_UNIVERSAL * T / (v - b) - a_alpha / (v**2 + 2 * b * v - b**2))[()]


def _roots(T, p, Tc, pc, omega):
    """Real compressibility roots above :math:`B` at ``(T, p)``, ascending, with ``A`` and ``B``."""
    a_alpha, b, _ = _constants(T, Tc, pc, omega)
    A, B = a_alpha * p / (R_UNIVERSAL * T) ** 2, b * p / (R_UNIVERSAL * T)
    z = np.roots([1.0, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)])
    real = np.sort(z[(np.abs(z.imag) < 1e-9 * np.abs(z.real).max()) & (z.real > B)].real)
    return real, A, B


def _log_fugacity_coefficient(Z, A, B):
    return Z - 1 - math.log(Z - B) - A / (2 * _SQRT2 * B) * math.log((Z + (1 + _SQRT2) * B) / (Z + (1 - _SQRT2) * B))


def peng_robinson_compressibility(T, p, Tc, pc, omega, phase="vapour"):
    r"""Compressibility factor :math:`Z = pv/(RT)` from the Peng-Robinson cubic.

    Parameters
    ----------
    T : float
        Absolute temperature [K], positive.
    p : float
        Pressure [Pa], positive.
    Tc, pc, omega : float
        Critical temperature [K], critical pressure [Pa] and acentric
        factor [-].
    phase : {"vapour", "liquid"}, optional
        Which root to return where the cubic has three: the largest
        (vapour) or the smallest (liquid). With one real root it is
        returned for either.

    Returns
    -------
    float
        Compressibility factor [-].

    Notes
    -----
    The cubic :math:`Z^3 - (1 - B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0`
    with :math:`A = a\alpha p/(RT)^2` and :math:`B = bp/(RT)` is solved with
    :func:`numpy.roots`. Which root is the stable phase is decided by the
    saturation pressure (:func:`peng_robinson_saturation_pressure`). Every
    parameter is scalar-only: a root is chosen for one state. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    D.-Y. Peng and D. B. Robinson (1976), as above.

    Examples
    --------
    >>> round(peng_robinson_compressibility(300.0, 101325.0, 190.56, 4.599e6, 0.011), 4)
    0.9978
    """
    _check(T, Tc, pc, omega)
    _positive(p=p)
    if np.ndim(T) or np.ndim(p):
        raise ValueError("T and p must be scalars")
    if phase not in ("vapour", "liquid"):
        raise ValueError("phase must be 'vapour' or 'liquid'")
    roots, _, _ = _roots(float(T), float(p), Tc, pc, omega)
    return float(roots[-1] if phase == "vapour" else roots[0])


def peng_robinson_saturation_pressure(T, Tc, pc, omega):
    r"""Vapour pressure from the Peng-Robinson equation: equal liquid and vapour fugacities.

    Parameters
    ----------
    T : float
        Absolute temperature [K], from :math:`0.45\,T_c` to
        :math:`0.9999\,T_c`.
    Tc, pc, omega : float
        Critical temperature [K], critical pressure [Pa] and acentric
        factor [-].

    Returns
    -------
    float
        Saturation pressure [Pa].

    Notes
    -----
    Between the two spinodal pressures, where :math:`(\partial p/\partial v)_T = 0`,
    the cubic has a liquid and a vapour root; Brent's method finds the
    pressure where their fugacity coefficients are equal. That is Maxwell's
    equal-area rule in another form. Every parameter is scalar-only: the
    search runs on one isotherm. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    D.-Y. Peng and D. B. Robinson (1976), as above.

    Examples
    --------
    >>> round(peng_robinson_saturation_pressure(0.7 * 190.56, 190.56, 4.599e6, 0.011) / 4.599e6, 4)
    0.0981
    """
    _check(T, Tc, pc, omega)
    if np.ndim(T):
        raise ValueError("T must be a scalar")
    _interval("T", T, _T_LOW * Tc, _T_TOP * Tc)
    return _saturation_pressure(float(T), float(Tc), float(pc), float(omega))


@lru_cache(maxsize=4096)
def _saturation_pressure(T, Tc, pc, omega):
    a_alpha, b, _ = _constants(T, Tc, pc, omega)
    RT = R_UNIVERSAL * T
    # spinodals: RT (v^2 + 2bv - b^2)^2 = 2 a alpha (v + b)(v - b)^2
    quad = np.polynomial.Polynomial([-(b**2), 2 * b, 1.0])
    poly = RT * quad**2 - 2 * a_alpha * np.polynomial.Polynomial([b, 1.0]) * np.polynomial.Polynomial([-b, 1.0]) ** 2
    v = np.sort(np.real([r for r in poly.roots() if abs(r.imag) < 1e-9 * abs(r.real) and r.real > b]))
    p_spin = RT / (v - b) - a_alpha / (v**2 + 2 * b * v - b**2)
    low, high = max(float(p_spin[0]), 1e-12 * pc), float(p_spin[1])

    def residual(p):
        roots, A, B = _roots(T, p, Tc, pc, omega)
        return _log_fugacity_coefficient(roots[0], A, B) - _log_fugacity_coefficient(roots[-1], A, B)

    # pull the brackets in by a hair so both roots still exist
    return brentq(residual, low * (1 + 1e-9), high * (1 - 1e-9), xtol=1e-12 * pc, rtol=1e-12)


def _saturation_temperature(fluid, p):
    T_low, T_top = _T_LOW * fluid.Tc, _T_TOP * fluid.Tc
    p_low, p_top = (_saturation_pressure(T, fluid.Tc, fluid.pc, fluid.omega) for T in (T_low, T_top))
    if not p_low <= p <= p_top:
        raise ValueError(f"saturation needs {p_low:.6g} Pa <= p <= {p_top:.6g} Pa (0.45 to 0.9999 Tc)")
    return brentq(lambda T: _saturation_pressure(T, fluid.Tc, fluid.pc, fluid.omega) - p, T_low, T_top, xtol=1e-10)


def _raw(fluid, T, p, phase):
    """Volume, enthalpy and entropy (ideal gas zero at 273.15 K and 101325 Pa) of one root."""
    roots, A, B = _roots(T, p, fluid.Tc, fluid.pc, fluid.omega)
    Z = roots[0] if phase == "liquid" else roots[-1]
    a_alpha, b, da = _constants(T, fluid.Tc, fluid.pc, fluid.omega)
    log_term = math.log((Z + (1 + _SQRT2) * B) / (Z + (1 - _SQRT2) * B))
    h_dep = R_UNIVERSAL * T * (Z - 1) + (T * da - a_alpha) / (2 * _SQRT2 * b) * log_term
    s_dep = R_UNIVERSAL * math.log(Z - B) + da / (2 * _SQRT2 * b) * log_term
    M = fluid.molar_mass
    h = fluid.cp0 * (T - 273.15) + h_dep / M
    s = fluid.cp0 * math.log(T / 273.15) - R_UNIVERSAL / M * math.log(p / 101325.0) + s_dep / M
    return Z * R_UNIVERSAL * T / (p * M), h, s


@lru_cache(maxsize=64)
def _offset(fluid):
    """Shifts that put saturated liquid at 0 °C at h = 200 kJ/kg, s = 1 kJ/(kg K) (IIR)."""
    p0 = _saturation_pressure(273.15, fluid.Tc, fluid.pc, fluid.omega)
    _, h, s = _raw(fluid, 273.15, p0, "liquid")
    return 200e3 - h, 1e3 - s


def _single(fluid, T, p, phase):
    v, h, s = _raw(fluid, T, p, "liquid" if phase == "liquid" else "vapour")
    dh, ds = _offset(fluid)
    return RefrigerantState(p=float(p), T=float(T), v=float(v), h=float(h + dh), s=float(s + ds), x=math.nan, phase=phase)


def _phase(fluid, T, p):
    if T >= fluid.Tc:
        return "supercritical"
    return "liquid" if p >= _saturation_pressure(T, fluid.Tc, fluid.pc, fluid.omega) else "vapour"


def _wet(fluid, p, x):
    T = _saturation_temperature(fluid, p)
    f, g = _single(fluid, T, p, "liquid"), _single(fluid, T, p, "vapour")
    mix = {k: (1 - x) * getattr(f, k) + x * getattr(g, k) for k in ("v", "h", "s")}
    return RefrigerantState(p=float(p), T=T, x=float(x), phase="mixture", **mix)


def refrigerant_properties(fluid, p, T=None, *, x=None, s=None, h=None):
    r"""State of a Peng-Robinson fluid from pressure and one other property.

    Parameters
    ----------
    fluid : CubicFluid
        The fluid, e.g. :data:`R134A`.
    p : float
        Pressure [Pa], positive.
    T : float, optional
        Temperature [K], from :math:`0.45\,T_c` to :math:`2.5\,T_c`, for a
        single-phase state.
    x : float, optional
        Vapour quality [-] in [0, 1] for a saturated mixture; needs ``p``
        on the saturation line between :math:`0.45\,T_c` and
        :math:`0.9999\,T_c` (up to about :math:`0.999\,p_c`).
    s : float, optional
        Specific entropy [J/(kg K)], e.g. after an isentropic compressor.
    h : float, optional
        Specific enthalpy [J/kg], e.g. after an expansion valve.

    Exactly one of ``T``, ``x``, ``s`` and ``h`` is given.

    Returns
    -------
    RefrigerantState
        Pressure, temperature, specific volume, enthalpy, entropy, quality
        (NaN outside the dome) and phase.

    Notes
    -----
    Enthalpy and entropy are the ideal gas with constant ``fluid.cp0``
    plus the Peng-Robinson departure functions

    .. math:: h - h^{ig} = RT(Z - 1) + \frac{T\,d(a\alpha)/dT - a\alpha}{2\sqrt2\,b}
              \ln\frac{Z + (1 + \sqrt2)B}{Z + (1 - \sqrt2)B},

    .. math:: s - s^{ig} = R\ln(Z - B) + \frac{d(a\alpha)/dT}{2\sqrt2\,b}
              \ln\frac{Z + (1 + \sqrt2)B}{Z + (1 - \sqrt2)B},

    shifted to the IIR reference (200 kJ/kg and 1 kJ/(kg K) for saturated
    liquid at 0 °C). Given ``T``, the phase is liquid at or above the
    saturation pressure and vapour below it. Given ``s`` or ``h``, the
    temperature is found by Brent's method. Peng-Robinson liquid volumes
    are typically 5-15% off; enthalpy differences are good to a few
    percent. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    D.-Y. Peng and D. B. Robinson (1976), as above.

    Examples
    --------
    >>> st = refrigerant_properties(R134A, 0.14e6, x=1.0)
    >>> round(st.T - 273.15, 1), st.phase
    (-18.7, 'mixture')
    """
    if not isinstance(fluid, CubicFluid):
        raise ValueError("fluid must be a CubicFluid")
    _positive(p=p)
    given = {name: value for name, value in (("T", T), ("x", x), ("s", s), ("h", h)) if value is not None}
    if len(given) != 1:
        raise ValueError("give exactly one of T, x, s or h")
    ((name, value),) = given.items()
    _finite(**{name: value})
    if np.ndim(p) or np.ndim(value):
        raise ValueError("p and the given property must be scalars")
    p, value = float(p), float(value)
    T_low, T_high = _T_LOW * fluid.Tc, _T_HIGH * fluid.Tc
    if name == "T":
        _interval("T", value, T_low, T_high)
        return _single(fluid, value, p, _phase(fluid, value, p))
    if name == "x":
        _interval("x", value, 0, 1)
        return _wet(fluid, p, value)
    if p < fluid.pc:
        T_sat = _saturation_temperature(fluid, p)
        f, g = (getattr(_single(fluid, T_sat, p, phase), name) for phase in ("liquid", "vapour"))
        if f <= value <= g:
            return _wet(fluid, p, (value - f) / (g - f))
        low, high, phase = (T_low, T_sat, "liquid") if value < f else (T_sat, T_high, "vapour")
    else:
        low, high, phase = T_low, T_high, "supercritical"
    edge_low, edge_high = (getattr(_single(fluid, T, p, phase), name) - value for T in (low, high))
    if edge_low > 0 or edge_high < 0:
        raise ValueError(f"{name} is outside the range covered at this pressure")
    T = brentq(lambda t: getattr(_single(fluid, t, p, phase), name) - value, low, high, xtol=1e-10)
    return _single(fluid, T, p, phase if phase != "supercritical" else _phase(fluid, T, p))
