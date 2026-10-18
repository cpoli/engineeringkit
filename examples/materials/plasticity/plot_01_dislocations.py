r"""
Taylor, Orowan and Polanyi: dislocations (1934)
===============================================

Jacob Frenkel (1926) estimated the shear stress needed to slide one plane of
a perfect crystal over the next as about :math:`G/2\pi`
(:func:`~engineeringkit.materials.ideal_shear_strength`). Single crystals of
copper and aluminium slip at stresses of order 1 MPa, thousands of times
less. In 1934 Geoffrey Taylor, Egon Orowan and Michael Polanyi independently
explained the gap with the dislocation, a line defect whose glide moves one
row of atoms at a time. The lattice resists it only with the Peierls stress
(Peierls 1940, Nabarro 1947), which falls exponentially as the slip planes
get farther apart (:func:`~engineeringkit.materials.peierls_stress`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import ideal_shear_strength, peierls_stress

# %%
# The strength gap
# ----------------
#
# Typical critical resolved shear stresses of annealed, high-purity single
# crystals against Frenkel's estimate.
# In a face-centred cubic metal the close-packed :math:`\{111\}` planes are
# :math:`a_0/\sqrt3` apart and the Burgers vector is :math:`a_0/\sqrt2`.

metals = {"copper": (48e9, 0.34, 0.361e-9, 0.65e6), "aluminium": (26e9, 0.35, 0.405e-9, 0.79e6)}
for name, (G, nu, a0, crss) in metals.items():
    b, d = a0 / np.sqrt(2), a0 / np.sqrt(3)
    ideal = ideal_shear_strength(G, b, d)
    tau_p = peierls_stress(G, nu, b, d)
    print(f"{name:9s}: ideal {ideal / 1e9:5.2f} GPa, Peierls {tau_p / 1e6:5.1f} MPa, measured {crss / 1e6:.2f} MPa")
    print(f"           ideal / measured = {ideal / crss:,.0f}")
    assert ideal / crss > 5000 and tau_p < ideal / 100

# %%
# Wide planes slip easily
# -----------------------
#
# The Peierls-Nabarro estimate is crude, it still exceeds the measured
# stress in FCC metals, whose dislocations split into partials and spread
# further, but its trend is right: closely spaced planes with short Burgers
# vectors resist slip far more than widely spaced close-packed ones. That is
# why slip happens on close-packed planes in close-packed directions, and
# why covalent crystals with narrow dislocations are brittle.

G, nu = 48e9, 0.34
ratio = np.linspace(0.5, 1.5, 101)
tau = peierls_stress(G, nu, 1.0, ratio)
assert np.all(np.diff(tau) < 0) and np.isclose(tau[0] / tau[-1], np.exp(2 * np.pi / (1 - nu)))

fig, ax = plt.subplots(figsize=(6, 3.8))
ax.semilogy(ratio, tau / G, label="Peierls-Nabarro")
ax.axhline(1 / (2 * np.pi), color="k", ls="--", label="Frenkel, G / 2 pi")
ax.axhline(0.65e6 / 48e9, color="C2", ls=":", label="copper, measured")
ax.axvline(np.sqrt(2 / 3), color="0.6", lw=0.8, label="FCC {111}")
ax.set(xlabel="plane spacing / Burgers vector", ylabel="shear stress / G", title="Why dislocations make metals soft")
ax.legend(loc="lower left")
fig.tight_layout()

plt.show()
