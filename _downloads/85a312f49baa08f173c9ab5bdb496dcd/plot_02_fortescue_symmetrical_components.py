r"""
Fortescue's symmetrical components (1918)
=========================================

An unbalanced three-phase system, such as one with a fault to ground,
couples all three phases and resisted hand calculation. Charles LeGeyt
Fortescue of Westinghouse showed in 1918 that any three phasors are the
sum of three *balanced* sets: a zero-sequence set (three equal phasors), a
positive-sequence set (a-b-c rotation) and a negative-sequence set (a-c-b)
(:func:`~engineeringkit.electrical_engineering.symmetrical_components`).
In a symmetrical network the three sequences do not interact, so each is
solved as a single-phase circuit and the results recombined
(:func:`~engineeringkit.electrical_engineering.phase_components`). This
example decomposes an unbalanced set, computes a single-line-to-ground
fault the way protection engineers still do and compares the four fault
types (:func:`~engineeringkit.electrical_engineering.fault_currents`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import fault_currents, phase_components, phasor, symmetrical_components
from engineeringkit.electrical_engineering.visualizers import plot_phasors

# %%
# Decomposing an unbalanced set
# -----------------------------

va, vb, vc = phasor(1.0, 0.0), phasor(0.8, np.radians(-125)), phasor(1.1, np.radians(115))
v0, v1, v2 = symmetrical_components(va, vb, vc)
print(f"|V0| = {abs(v0):.3f}, |V1| = {abs(v1):.3f}, |V2| = {abs(v2):.3f} pu; unbalance V2/V1 = {abs(v2 / v1):.1%}")
zero, positive, negative = phase_components(v0, 0, 0), phase_components(0, v1, 0), phase_components(0, 0, v2)
assert np.allclose(np.add(np.add(zero, positive), negative), (va, vb, vc))

# %%
# A single-line-to-ground fault
# -----------------------------
#
# Phase a shorts to ground at a bus whose sequence impedances are
# :math:`Z_1 = Z_2 = j0.2` and :math:`Z_0 = j0.1` pu, prefault voltage 1 pu.
# The boundary conditions :math:`I_b = I_c = 0`, :math:`V_a = 0` connect the
# three sequence networks in series: :math:`I_0 = I_1 = I_2 = E/(Z_0 + Z_1 + Z_2)`.

E, Z1, Z2, Z0 = 1.0, 0.2j, 0.2j, 0.1j
I_seq = E / (Z0 + Z1 + Z2)
Ia, Ib, Ic = phase_components(I_seq, I_seq, I_seq)
V_seq = (-Z0 * I_seq, E - Z1 * I_seq, -Z2 * I_seq)
Va, Vb, Vc = phase_components(*V_seq)
print(f"fault current |Ia| = {abs(Ia):.2f} pu (a three-phase fault would give {abs(E / Z1):.2f} pu)")
print(f"healthy-phase voltages |Vb| = {abs(Vb):.3f}, |Vc| = {abs(Vc):.3f} pu")
assert np.isclose(Ia, 3 * E / (Z0 + Z1 + Z2))
assert abs(Ib) < 1e-12 and abs(Ic) < 1e-12 and abs(Va) < 1e-12
assert np.allclose(fault_currents(E, Z1, Z2, Z0, "slg"), (Ia, Ib, Ic))

# %%
# The four fault types
# --------------------
#
# The other faults connect the sequence networks differently. Here the
# single-line-to-ground fault is the most severe, because the zero-sequence
# impedance is below the positive-sequence one.

for kind in ("3ph", "slg", "ll", "llg"):
    currents = fault_currents(E, Z1, Z2, Z0, kind)
    print(f"{kind:>4}: |Ia|, |Ib|, |Ic| = " + ", ".join(f"{abs(i):.2f}" for i in currents) + " pu")
assert np.isclose(abs(fault_currents(E, Z1, Z2, Z0, "ll")[1]), np.sqrt(3) / 2 * abs(E / Z1))

fig, axes = plt.subplots(1, 4, figsize=(14, 3.6))
labels = ["a", "b", "c"]
plot_phasors([va, vb, vc], labels=labels, ax=axes[0])
plot_phasors(positive, labels=labels, ax=axes[1])
plot_phasors(negative, labels=labels, ax=axes[2])
plot_phasors(zero, labels=labels, ax=axes[3])
for ax, title in zip(axes, ["unbalanced set", "positive sequence", "negative sequence", "zero sequence"]):
    ax.set_title(title)
    ax.set(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
fig.tight_layout()

plt.show()
