r"""
Dalton's law of partial pressures (1801)
========================================

John Dalton, a Manchester schoolteacher, read a set of essays to the city's
Literary and Philosophical Society in October 1801. He argued that each gas
in a mixture presses on the walls as if the others were absent, so the
total pressure of a mixture is the sum of the partial pressures,
:math:`p = \sum p_i` with :math:`p_i = y_i p`
(:func:`~engineeringkit.thermodynamics.partial_pressures`). The rule underlies
moist-air calculations, combustion products and every ideal-gas mixture: its
molar mass is the mole-weighted mean
(:func:`~engineeringkit.thermodynamics.mixture_molar_mass`), and mixing
different gases generates entropy
(:func:`~engineeringkit.thermodynamics.entropy_of_mixing`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import (
    R_UNIVERSAL,
    entropy_change_isothermal_ideal_gas,
    entropy_of_mixing,
    mass_fractions,
    mixture_molar_mass,
    mole_fractions,
    partial_pressures,
)

# %%
# Dry air from its composition
# ----------------------------
#
# By volume (and so by moles) dry air is 78.08% nitrogen, 20.95% oxygen and
# 0.93% argon, plus traces. Dalton's law gives each gas its share of the
# atmosphere's pressure; the oxygen partial pressure is what a climber's lungs
# feel.

names = ["N2", "O2", "Ar"]
M = np.array([0.028013, 0.031999, 0.039948])
y = np.array([0.7808, 0.2095, 0.0097])
y = y / y.sum()
p_i = partial_pressures(y, 101325.0)
x = mass_fractions(y, M)
for name, yi, xi, pi in zip(names, y, x, p_i):
    print(f"{name}: mole fraction {yi:.4f}, mass fraction {xi:.4f}, partial pressure {pi / 1e3:.2f} kPa")
M_air = mixture_molar_mass(y, M)
print(f"apparent molar mass {M_air * 1e3:.3f} g/mol, gas constant {R_UNIVERSAL / M_air:.1f} J/(kg K)")
assert np.isclose(p_i.sum(), 101325.0) and np.isclose(M_air, 0.02897, rtol=1e-3)
assert np.allclose(mole_fractions(x, M), y)

# %%
# The entropy of mixing
# ---------------------
#
# Remove a partition between two different gases at the same temperature and
# pressure: each expands into the whole volume, down to its partial pressure.
# The entropy rises by :math:`-nR\sum y_i\ln y_i`, the sum of the two
# isothermal expansions, and is largest for equal amounts.

n = 1.0
by_expansion = sum(entropy_change_isothermal_ideal_gas(n * yi, 300.0, yi, 1.0) for yi in y)
print(f"mixing 1 mol of air from its pure components: {entropy_of_mixing(y, n):.3f} J/K (expansions {by_expansion:.3f} J/K)")
assert np.isclose(entropy_of_mixing(y, n), by_expansion)

share = np.linspace(0, 1, 101)
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(share, [entropy_of_mixing([s, 1 - s]) for s in share])
ax.axhline(R_UNIVERSAL * np.log(2), color="gray", ls=":", label="R ln 2")
ax.set(xlabel="mole fraction of the first gas", ylabel="entropy of mixing [J/(mol K)]", title="Mixing two ideal gases")
ax.legend()
fig.tight_layout()
assert np.isclose(max(entropy_of_mixing([s, 1 - s]) for s in share), R_UNIVERSAL * np.log(2))

plt.show()
