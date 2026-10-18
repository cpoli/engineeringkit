r"""
Thévenin's equivalent circuit
=============================

Léon Charles Thévenin, an engineer of the French telegraph service, showed
in 1883 that any linear network seen from two terminals behaves like a
single EMF in series with a single resistance (Hermann von Helmholtz had
reached the same result in 1853). The EMF is the open-circuit voltage; the
resistance is that voltage divided by the short-circuit current, or the
resistance seen with every source zeroed. This example reduces a two-source
network with :func:`~engineeringkit.electrical_engineering.nodal_analysis`,
builds the equivalent with
:func:`~engineeringkit.electrical_engineering.thevenin_from_source`, and
checks that the two agree for every load.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import nodal_analysis, thevenin_from_source

# %%
# The network and its terminals
# -----------------------------
#
# A 10 V source and a 5 mA source feed a resistive network; the load is
# connected between node 3 and ground.

network = [(1, 2, 2e3), (2, 0, 3e3), (2, 3, 1e3)]
vsrc, isrc = [(1, 0, 10.0)], [(0, 3, 5e-3)]


def solve(load=None):
    resistors = network + ([] if load is None else [(3, 0, load)])
    return nodal_analysis(3, resistors, current_sources=isrc, voltage_sources=vsrc)


v_open = solve().voltages[3]
i_short = solve(1e-6).voltages[3] / 1e-6  # a 1 µohm load is a short
vth, rth = thevenin_from_source(v_open, v_open / i_short)
# The same resistance with both sources zeroed: 1 k + (2 k || 3 k).
assert np.isclose(rth, 1e3 + 2e3 * 3e3 / 5e3, rtol=1e-6)
print(f"V_th = {vth:.3f} V, R_th = {rth:.1f} ohm")

# %%
# Any load sees the same thing
# ----------------------------

loads = np.logspace(1, 5, 25)
v_network = np.array([solve(r).voltages[3] for r in loads])
v_thevenin = vth * loads / (rth + loads)
assert np.allclose(v_network, v_thevenin, rtol=1e-6)

fig, ax = plt.subplots()
ax.semilogx(loads, v_network, "o", label="full network")
ax.semilogx(loads, v_thevenin, label="Thévenin equivalent")
ax.set(xlabel="load resistance [ohm]", ylabel="terminal voltage [V]", title="Thévenin's theorem")
ax.legend()

plt.show()
