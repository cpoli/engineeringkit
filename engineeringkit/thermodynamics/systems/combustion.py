r"""Combustion: stoichiometry, heating values and the adiabatic flame temperature.

A hydrocarbon fuel :math:`\mathrm{C}_x\mathrm{H}_y\mathrm{O}_z` burns completely
in air (21% oxygen and 79% nitrogen by volume, 3.76 mol of nitrogen per mol of
oxygen) to carbon dioxide and water,

.. math:: \mathrm{C}_x\mathrm{H}_y\mathrm{O}_z + (1 + e)\,a\,(\mathrm{O}_2 + 3.76\,\mathrm{N}_2)
   \to x\,\mathrm{CO}_2 + \tfrac{y}{2}\mathrm{H}_2\mathrm{O} + e\,a\,\mathrm{O}_2
   + 3.76(1 + e)\,a\,\mathrm{N}_2,

with :math:`a = x + y/4 - z/2` and excess air :math:`e`. Hess's law makes the
heat of reaction the difference of the enthalpies of formation of products and
reactants. When none of that heat leaves, the products reach the adiabatic
flame temperature. Product enthalpies come from NIST's Shomate equations
(:mod:`~engineeringkit.thermodynamics.utils.shomate`). Dissociation is
neglected, so the flame temperature is an upper bound.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.thermodynamics.utils.shomate import FORMATION_ENTHALPY, LIQUID_WATER_FORMATION_ENTHALPY, sensible_enthalpy

__all__ = [
    "MOLAR_MASS_AIR",
    "fuel_molar_mass",
    "stoichiometric_air_fuel_ratio",
    "combustion_products",
    "heating_value",
    "adiabatic_flame_temperature",
]

MOLAR_MASS_AIR = 0.02897
"""Molar mass of dry air [kg/mol]."""

_N2_PER_O2 = 3.76
_ATOMIC_MASS = {"C": 0.012011, "H": 0.001008, "O": 0.015999}  # kg/mol


def _oxygen(c, h, o):
    """Validate the formula and return the theoretical oxygen :math:`x + y/4 - z/2` [mol/mol fuel]."""
    _nonnegative(c=c, h=h, o=o)
    need = c + h / 4 - o / 2
    _positive(**{"c + h/4 - o/2 (oxygen needed)": need})
    return need


def fuel_molar_mass(c, h, o=0.0):
    r"""Molar mass of a fuel :math:`\mathrm{C}_c\mathrm{H}_h\mathrm{O}_o`.

    Parameters
    ----------
    c, h, o : float
        Atoms of carbon, hydrogen and oxygen per molecule [-], nonnegative.

    Returns
    -------
    float
        Molar mass [kg/mol].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(fuel_molar_mass(1, 4), 5)
    0.01604
    """
    _nonnegative(c=c, h=h, o=o)
    return c * _ATOMIC_MASS["C"] + h * _ATOMIC_MASS["H"] + o * _ATOMIC_MASS["O"]


def stoichiometric_air_fuel_ratio(c, h, o=0.0):
    r"""Mass of air per mass of fuel for complete combustion with no excess air.

    Parameters
    ----------
    c, h, o : float
        Atoms of carbon, hydrogen and oxygen per fuel molecule [-],
        nonnegative, needing some oxygen.

    Returns
    -------
    float
        :math:`4.76\,a\,M_{air}/M_{fuel}` [kg air/kg fuel], with
        :math:`a = c + h/4 - o/2` mol of oxygen per mol of fuel.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(stoichiometric_air_fuel_ratio(8, 18), 2)  # octane
    15.09
    """
    need = _oxygen(c, h, o)
    return need * (1 + _N2_PER_O2) * MOLAR_MASS_AIR / fuel_molar_mass(c, h, o)


def combustion_products(c, h, o=0.0, excess_air=0.0):
    r"""Moles of each product of complete combustion per mole of fuel.

    Parameters
    ----------
    c, h, o : float or array_like
        Atoms of carbon, hydrogen and oxygen per fuel molecule [-].
    excess_air : float or array_like, optional
        Air beyond the stoichiometric amount, as a fraction [-],
        nonnegative (0.5 is 150% theoretical air).

    Returns
    -------
    dict
        Moles of ``"CO2"``, ``"H2O"``, ``"O2"`` and ``"N2"`` [mol/mol fuel].

    Notes
    -----
    Rich mixtures, which also produce CO and H2, are not covered. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> {k: round(float(v), 2) for k, v in combustion_products(1, 4).items()}
    {'CO2': 1.0, 'H2O': 2.0, 'O2': 0.0, 'N2': 7.52}
    """
    need = _oxygen(c, h, o)
    _nonnegative(excess_air=excess_air)
    c, h, excess_air, need = (np.asarray(v, dtype=float) for v in (c, h, excess_air, need))
    return {
        "CO2": (c + 0 * need)[()],
        "H2O": (h / 2 + 0 * need)[()],
        "O2": (excess_air * need)[()],
        "N2": (_N2_PER_O2 * (1 + excess_air) * need)[()],
    }


def heating_value(c, h, o=0.0, *, h_f_fuel, water="vapour"):
    r"""Heating value of a fuel from enthalpies of formation (Hess's law).

    Parameters
    ----------
    c, h, o : float
        Atoms of carbon, hydrogen and oxygen per fuel molecule [-].
    h_f_fuel : float
        Standard enthalpy of formation of the fuel at 298.15 K [J/mol], in
        the phase in which it is burned (e.g. -74.87e3 for methane gas,
        -249.95e3 for liquid octane).
    water : {"vapour", "liquid"}, optional
        Phase of the product water: vapour for the lower heating value,
        liquid for the higher.

    Returns
    -------
    float
        Heat released per kilogram of fuel burned at 25 °C with its products
        cooled back to 25 °C [J/kg].

    Notes
    -----
    :math:`HV = [h_f^{fuel} - c\,h_f^{CO_2} - (h/2)\,h_f^{H_2O}]/M_{fuel}`;
    the excess air drops out because oxygen and nitrogen have zero enthalpy
    of formation. Inputs are raw SI numbers with no dimensional metadata;
    see :doc:`/guides/thermodynamics`.

    References
    ----------
    G. H. Hess, "Thermochemische Untersuchungen," *Ann. Phys.* 126 (1840)
    385-404, doi:10.1002/andp.18401260620.

    Examples
    --------
    >>> round(heating_value(1, 4, h_f_fuel=-74.87e3) / 1e6, 2)
    50.01
    """
    _oxygen(c, h, o)
    _finite(h_f_fuel=h_f_fuel)
    if water not in ("vapour", "liquid"):
        raise ValueError("water must be 'vapour' or 'liquid'")
    h_f_water = FORMATION_ENTHALPY["H2O"] if water == "vapour" else LIQUID_WATER_FORMATION_ENTHALPY
    released = h_f_fuel - c * FORMATION_ENTHALPY["CO2"] - h / 2 * h_f_water
    return released / fuel_molar_mass(c, h, o)


def adiabatic_flame_temperature(c, h, o=0.0, *, h_f_fuel, excess_air=0.0, T_air=298.15):
    r"""Temperature of the products of complete, adiabatic, constant-pressure combustion.

    Parameters
    ----------
    c, h, o : float
        Atoms of carbon, hydrogen and oxygen per fuel molecule [-].
    h_f_fuel : float
        Standard enthalpy of formation of the fuel at 298.15 K [J/mol].
    excess_air : float, optional
        Excess air as a fraction of the stoichiometric amount [-],
        nonnegative.
    T_air : float, optional
        Temperature of the combustion air [K], 298.15 K to 1500 K; the fuel
        enters at 298.15 K.

    Returns
    -------
    float
        Adiabatic flame temperature [K].

    Notes
    -----
    Solves :math:`\sum_p n_i[h_{f,i} + \Delta h_i(T)] = h_f^{fuel} + \sum_r n_j\Delta h_j(T_{air})`
    by Brent's method, with water as vapour and no dissociation. Above
    about 2000 K real flames dissociate CO2 and H2O, so the measured
    temperature is lower (methane in air: about 2230 K at equilibrium
    against 2327 K here). Every parameter is scalar-only: the balance is
    solved for one fuel and mixture. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    Y. A. Çengel and M. A. Boles, *Thermodynamics: An Engineering
    Approach*, 9th ed., New York: McGraw-Hill (2019), Ch. 15.

    Examples
    --------
    >>> round(adiabatic_flame_temperature(1, 4, h_f_fuel=-74.87e3))
    2327
    """
    products = combustion_products(c, h, o, excess_air)
    _finite(h_f_fuel=h_f_fuel)
    _interval("T_air", T_air, 298.15, 1500.0)
    need = c + h / 4 - o / 2
    reactants = h_f_fuel + (1 + excess_air) * need * (sensible_enthalpy("O2", T_air) + _N2_PER_O2 * sensible_enthalpy("N2", T_air))
    formation = sum(n * FORMATION_ENTHALPY[k] for k, n in products.items())

    def residual(T):
        return formation + sum(n * sensible_enthalpy(k, T) for k, n in products.items()) - reactants

    if residual(298.15) > 0:
        raise ValueError("the reaction does not release heat")
    return float(brentq(residual, 298.15, 6000.0, xtol=1e-8))
