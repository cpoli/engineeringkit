r"""
Wheatstone's bridge (1843)
==========================

Samuel Hunter Christie described the circuit in 1833, and Charles
Wheatstone made it the standard way to measure resistance in his 1843
Bakerian lecture. Two voltage dividers share a source; a galvanometer joins
their midpoints. When :math:`R_1/R_2 = R_3/R_4` the midpoints sit at the same
potential and the galvanometer reads zero, whatever the source voltage
(:func:`~engineeringkit.electrical_engineering.wheatstone_balanced`). A null
reading needs no calibrated meter, which is why the bridge is so accurate.
This example balances a bridge against an unknown resistor and shows how
sharply the null depends on it.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import nodal_analysis, wheatstone_balanced

# %%
# The galvanometer current
# ------------------------
#
# Nodes: 1 the source, 2 and 3 the midpoints. The unknown :math:`R_x` is
# :math:`R_4`; :math:`R_3` is the calibrated adjustable arm.

R1, R2, Rx, Rg, V = 1000.0, 1000.0, 472.0, 50.0, 5.0


def galvanometer_current(R3):
    bridge = [(1, 2, R1), (2, 0, R2), (1, 3, R3), (3, 0, Rx), (2, 3, Rg)]
    return nodal_analysis(3, bridge, voltage_sources=[(1, 0, V)]).resistor_currents[4]


R3 = np.linspace(400, 550, 301)
ig = np.array([galvanometer_current(r) for r in R3])
balance = R3[np.argmin(np.abs(ig))]
print(f"null at R3 = {balance:.1f} ohm, so Rx = R3 R2 / R1 = {balance * R2 / R1:.1f} ohm")
assert wheatstone_balanced(R1, R2, balance, Rx)

# %%
# The null does not depend on the source
# --------------------------------------

for v_source in (1.0, 5.0, 20.0):
    bridge = [(1, 2, R1), (2, 0, R2), (1, 3, Rx), (3, 0, Rx), (2, 3, Rg)]
    i_null = nodal_analysis(3, bridge, voltage_sources=[(1, 0, v_source)]).resistor_currents[4]
    assert abs(i_null) < 1e-15

fig, ax = plt.subplots()
ax.plot(R3, ig * 1e6)
ax.axhline(0, color="gray", lw=0.5)
ax.axvline(Rx, color="C3", ls="--", label="$R_3 = R_x R_1/R_2$")
ax.set(xlabel="adjustable arm $R_3$ [ohm]", ylabel="galvanometer current [µA]", title="Balancing a Wheatstone bridge")
ax.legend()

plt.show()
