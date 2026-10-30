r"""
Bakhmeteff's specific energy (1912)
===================================

Boris Bakhmeteff, in a 1912 book on varied flow written at St Petersburg,
measured the energy of a stream from the channel bed instead of from a
fixed datum: the specific energy :math:`E = y + q^2/(2gy^2)`
(:func:`~engineeringkit.fluid_engineering.specific_energy`). Plotted
against depth, it has a least value at the critical depth
:math:`y_c = (q^2/g)^{1/3}`
(:func:`~engineeringkit.fluid_engineering.critical_depth`), where the
Froude number is 1, and any larger energy can be carried at a subcritical
and a supercritical alternate depth
(:func:`~engineeringkit.fluid_engineering.alternate_depths`). This example
draws the specific-energy diagram and uses it to find the flow over a
raised step in the bed.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import alternate_depths, critical_depth, froude_number, specific_energy

# %%
# The specific-energy diagram
# ---------------------------

y = np.linspace(0.08, 3.0, 300)
fig, ax = plt.subplots(figsize=(5, 5))
for q in [1.0, 2.0, 4.0]:
    yc = critical_depth(q)
    ax.plot(specific_energy(y, q), y, label=f"q = {q:g} m²/s")
    ax.plot(1.5 * yc, yc, "ko")
ax.plot([0, 3], [0, 3], "k:", lw=0.8, label="E = y")
ax.set(xlim=(0, 3), ylim=(0, 3), xlabel="specific energy E [m]", ylabel="depth y [m]", title="Specific energy")
ax.legend()
fig.tight_layout()

# %%
# Flow over a raised step
# -----------------------
#
# A stream 1.8 m deep carrying 2 m²/s meets a step 0.2 m high. With no
# loss, the specific energy over the step is :math:`E_1 - \Delta z`, and
# the subcritical stream stays subcritical: its surface dips.

q, y1, dz = 2.0, 1.8, 0.2
E1 = specific_energy(y1, q)
deep, shallow = alternate_depths(E1 - dz, q)
print(f"upstream: y = {y1} m, E = {E1:.3f} m, Fr = {froude_number(q / y1, y1):.2f}")
print(f"over the step: y = {deep:.3f} m (the supercritical alternate would be {shallow:.3f} m)")
print(f"water surface drops by {y1 - (deep + dz):.3f} m")

# %%
# The largest step the stream can pass without backing up is the one that
# brings the energy down to its minimum, :math:`3y_c/2`.

yc = critical_depth(q)
dz_max = E1 - 1.5 * yc
print(f"critical depth {yc:.3f} m; largest step {dz_max:.3f} m")

assert np.isclose(specific_energy(deep, q), E1 - dz) and deep > yc > shallow
assert np.isclose(froude_number(q / yc, yc), 1.0)
assert np.allclose(alternate_depths(E1 - dz_max, q), yc)

plt.show()
