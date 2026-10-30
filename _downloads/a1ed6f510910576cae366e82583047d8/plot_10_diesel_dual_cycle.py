r"""
Diesel's cycle and the dual cycle (1893)
========================================

Diesel's ideal engine burns its fuel at constant pressure. Real high-speed
diesels burn part of the charge almost at once, at nearly constant volume,
and the rest as the piston moves. The dual (limited-pressure) cycle models
this with both heat additions
(:func:`~engineeringkit.thermodynamics.dual_efficiency`). It contains Otto's
cycle (no constant-pressure part) and Diesel's (no constant-volume part) as
its limits.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import diesel_efficiency, dual_efficiency, isentropic_pressure_volume, otto_efficiency
from engineeringkit.thermodynamics.visualizers import plot_pv_diagram

# %%
# Four states of a dual cycle
# ---------------------------
#
# Compression ratio 16, pressure ratio 1.5 at constant volume, then cutoff
# ratio 1.4 at constant pressure, air with :math:`\gamma = 1.4`.

r, alpha, rc, g, T1, p1 = 16.0, 1.5, 1.4, 1.4, 300.0, 100e3
T2 = T1 * r ** (g - 1)
Tx, T3 = alpha * T2, alpha * T2 * rc
T4 = T3 * (rc / r) ** (g - 1)
eta_states = 1 - (T4 - T1) / ((Tx - T2) + g * (T3 - Tx))
print(f"T2 = {T2:.0f} K, Tx = {Tx:.0f} K, T3 = {T3:.0f} K, T4 = {T4:.0f} K")
print(f"efficiency from the states {eta_states:.4f}, from the formula {dual_efficiency(r, alpha, rc):.4f}")
assert np.isclose(eta_states, dual_efficiency(r, alpha, rc))
assert np.isclose(dual_efficiency(r, 1.0, rc), diesel_efficiency(r, rc)) and np.isclose(dual_efficiency(r, alpha, 1.0), otto_efficiency(r))

v1 = 1.0
v2 = v1 / r
v3 = rc * v2
compression = np.linspace(v1, v2, 40)
expansion = np.linspace(v3, v1, 40)
p2 = p1 * r**g
p3 = alpha * p2
v = np.r_[compression, v2, v3, expansion]
p = np.r_[isentropic_pressure_volume(p1, v1, compression, g), p3, p3, isentropic_pressure_volume(p3, v3, expansion, g)]
ax = plot_pv_diagram(np.r_[v, v1], np.r_[p, p1], label="dual cycle")
ax.set(xlabel="volume / initial volume", title="Dual cycle, r = 16")

# %%
# The same heat in different ways
# -------------------------------
#
# For a fixed compression ratio and heat input, adding more of the heat at
# constant volume raises the efficiency: Otto beats dual beats Diesel. Diesel
# engines win in practice because they can run at twice the compression ratio.

alphas = np.linspace(1.0, 2.5, 31)
fig, ax = plt.subplots(figsize=(6, 3.5))
heat = (alpha - 1) + g * alpha * (rc - 1)  # the example's heat input, in units of cv T2
cutoffs = 1 + (heat - (alphas - 1)) / (g * alphas)
valid = cutoffs >= 1
ax.plot(alphas[valid], dual_efficiency(r, alphas[valid], cutoffs[valid]), label="dual, same heat input")
ax.axhline(otto_efficiency(r), color="gray", ls=":", label="Otto")
ax.set(xlabel="constant-volume pressure ratio α", ylabel="efficiency", title="Splitting the heat addition, r = 16")
ax.legend()
fig.tight_layout()
assert np.all(np.diff(dual_efficiency(r, alphas[valid], cutoffs[valid])) > 0)

plt.show()
