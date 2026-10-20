r"""
Voigt and Reuss: bounds on a composite's stiffness
==================================================

How stiff is a mixture of two phases? Woldemar Voigt (1889) averaged the
elastic constants of the crystals in a polycrystal on the assumption that
all of them strain alike; Adolf Reuss (1929) averaged the compliances on the
assumption that all carry the same stress. For a fibre composite the two
assumptions are exact along and across the fibres
(:func:`~engineeringkit.materials.rule_of_mixtures_E_longitudinal`,
:func:`~engineeringkit.materials.rule_of_mixtures_E_transverse`), and Hill
(1952) proved that every other arrangement lies between them. The example
builds both from springs and checks a glass-epoxy laminate.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import density_rule_of_mixtures, rule_of_mixtures_E_longitudinal, rule_of_mixtures_E_transverse

# %%
# Springs in parallel and in series
# ---------------------------------
#
# A unit cube with a fibre fraction :math:`V_f` is a slab of fibre and a
# slab of matrix. Loaded along the slabs, they are springs in parallel with
# stiffnesses :math:`E_f V_f` and :math:`E_m(1 - V_f)`; loaded across them,
# springs in series with :math:`E_f/V_f` and :math:`E_m/(1 - V_f)`.

Ef, Em, Vf = 72e9, 3.5e9, 0.55
parallel = Ef * Vf + Em * (1 - Vf)
series = 1 / (1 / (Ef / Vf) + 1 / (Em / (1 - Vf)))
E1, E2 = rule_of_mixtures_E_longitudinal(Ef, Em, Vf), rule_of_mixtures_E_transverse(Ef, Em, Vf)
rho = density_rule_of_mixtures(2540.0, 1200.0, Vf)
print(f"glass-epoxy at Vf = {Vf}: E1 = {E1 / 1e9:.1f} GPa, E2 = {E2 / 1e9:.2f} GPa, density {rho:.0f} kg/m^3")
assert np.isclose(E1, parallel) and np.isclose(E2, series)

# %%
# The bounds and the gap between them
# -----------------------------------
#
# Any microstructure, particles, random fibres or foam, has a modulus in
# the shaded band. The band is wide when the phases differ a lot: at
# :math:`V_f = 0.5` the bounds differ by a factor of about 6.

fractions = np.linspace(0, 1, 201)
voigt = rule_of_mixtures_E_longitudinal(Ef, Em, fractions)
reuss = rule_of_mixtures_E_transverse(Ef, Em, fractions)
hill = (voigt + reuss) / 2
assert np.all(reuss <= hill) and np.all(hill <= voigt) and np.all(np.diff(voigt) > 0)
print(f"Voigt / Reuss at Vf = 0.5: {voigt[100] / reuss[100]:.1f}")

fig, ax = plt.subplots(figsize=(6, 3.8))
ax.fill_between(fractions, reuss / 1e9, voigt / 1e9, alpha=0.2, label="attainable")
ax.plot(fractions, voigt / 1e9, label="Voigt (equal strain)")
ax.plot(fractions, reuss / 1e9, label="Reuss (equal stress)")
ax.plot(fractions, hill / 1e9, "k--", lw=1, label="Hill average")
ax.set(xlabel="fibre volume fraction", ylabel="modulus [GPa]", title="Glass fibre in epoxy")
ax.legend()
fig.tight_layout()

plt.show()
