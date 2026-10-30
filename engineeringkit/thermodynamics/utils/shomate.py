r"""Ideal-gas enthalpies of the combustion gases from NIST's Shomate equations.

Each species has one coefficient set per temperature range; with
:math:`t = T/1000`,

.. math:: H^\circ - H^\circ_{298.15} = At + Bt^2/2 + Ct^3/3 + Dt^4/4 - E/t + F - H \quad [\mathrm{kJ/mol}].

Coefficients and standard enthalpies of formation are those of the NIST
Chemistry WebBook (Chase, *NIST-JANAF Thermochemical Tables*, 1998; CODATA
key values, Cox, Wagman and Medvedev, 1984).
"""

from __future__ import annotations

import numpy as np

__all__ = ["SPECIES", "FORMATION_ENTHALPY", "LIQUID_WATER_FORMATION_ENTHALPY", "sensible_enthalpy"]

# species: [(upper limit of the range [K], (A, B, C, D, E, F, G, H)), ...]; the lowest range extends down to 298.15 K
_SHOMATE = {
    "CO2": [
        (1200.0, (24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224)),
        (6000.0, (58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393.5224)),
    ],
    "H2O": [
        (1700.0, (30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, 223.3967, -241.8264)),
        (6000.0, (41.96426, 8.622053, -1.499780, 0.098119, -11.15764, -272.1797, 219.7809, -241.8264)),
    ],
    "N2": [
        (500.0, (28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 226.4168, 0.0)),
        (2000.0, (19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900, 0.0)),
        (6000.0, (35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 224.9810, 0.0)),
    ],
    "O2": [
        (700.0, (31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0)),
        (2000.0, (30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0)),
        (6000.0, (20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 237.6185, 0.0)),
    ],
}

SPECIES = tuple(_SHOMATE)
"""tuple of str: Species with Shomate data: CO2, H2O (gas), N2 and O2."""

FORMATION_ENTHALPY = {"CO2": -393.51e3, "H2O": -241.826e3, "N2": 0.0, "O2": 0.0}
"""dict: Standard enthalpy of formation at 298.15 K [J/mol], water as vapour (CODATA)."""

LIQUID_WATER_FORMATION_ENTHALPY = -285.83e3
"""float: Standard enthalpy of formation of liquid water at 298.15 K [J/mol] (CODATA)."""


def sensible_enthalpy(species, T):
    """Ideal-gas enthalpy above 298.15 K, :math:`H^\\circ(T) - H^\\circ(298.15)` [J/mol].

    Parameters
    ----------
    species : str
        One of :data:`SPECIES`.
    T : float or array_like
        Temperature [K], 298.15 K to 6000 K.

    Returns
    -------
    float or ndarray
        Sensible enthalpy [J/mol].
    """
    T = np.asarray(T, dtype=float)
    ranges = _SHOMATE[species]
    out = np.empty_like(T)
    lower = -np.inf
    for upper, (A, B, C, D, E, F, _, H) in ranges:
        mask = (T > lower) & (T <= upper)
        t = T[mask] / 1000
        out[mask] = 1e3 * (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F - H)
        lower = upper
    return out[()]
