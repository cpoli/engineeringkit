r"""
Norton's equivalent circuit
===========================

In 1926 Edward Lawry Norton at Bell Labs, and independently Hans Ferdinand
Mayer at Siemens, gave the dual of Thévenin's theorem: a linear two-terminal
network also behaves like a current source :math:`I_N = V_{th}/R_{th}` in
parallel with the same resistance
(:func:`~engineeringkit.electrical_engineering.norton_from_thevenin`).
Norton wrote it in an internal memorandum on filter design; Mayer used it
to model the vacuum-tube amplifier, which acts more like a current source
than a voltage source. This example checks that the Thévenin and Norton
forms are indistinguishable at the terminals, and that converting sources
back and forth simplifies a ladder to a single source.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import nodal_analysis, norton_from_thevenin, parallel_resistance, thevenin_from_source

# %%
# Same terminal behaviour
# -----------------------

vth, rth = thevenin_from_source(9.0, 300.0)
i_n, r_n = norton_from_thevenin(vth, rth)
print(f"Thévenin {vth:.1f} V + {rth:.0f} ohm  <->  Norton {i_n * 1e3:.1f} mA || {r_n:.0f} ohm")
loads = np.logspace(0, 4, 30)
v_thevenin = [nodal_analysis(2, [(1, 2, rth), (2, 0, r)], voltage_sources=[(1, 0, vth)]).voltages[2] for r in loads]
v_norton = [nodal_analysis(1, [(1, 0, r_n), (1, 0, r)], current_sources=[(0, 1, i_n)]).voltages[1] for r in loads]
assert np.allclose(v_thevenin, v_norton)

# %%
# Source transformations
# ----------------------
#
# Two batteries with internal resistance in parallel: convert each to
# Norton form, add the currents and the conductances, convert back.

batteries = [(12.0, 0.5), (11.4, 0.3)]
nortons = [norton_from_thevenin(v, r) for v, r in batteries]
i_total = sum(i for i, _ in nortons)
r_total = parallel_resistance(*(r for _, r in nortons))
v_eq = i_total * r_total
direct = nodal_analysis(3, [(1, 3, 0.5), (2, 3, 0.3)], voltage_sources=[(1, 0, 12.0), (2, 0, 11.4)]).voltages[3]
print(f"open-circuit voltage of the pair: {v_eq:.4f} V (nodal analysis {direct:.4f} V), R = {r_total:.4f} ohm")
assert np.isclose(v_eq, direct)

fig, ax = plt.subplots()
ax.semilogx(loads, v_thevenin, label="Thévenin form")
ax.semilogx(loads, v_norton, "--", label="Norton form")
ax.set(xlabel="load resistance [ohm]", ylabel="load voltage [V]", title="Norton's equivalent")
ax.legend()

plt.show()
