r"""
Brayton's constant-pressure cycle (1872)
========================================

George Brayton's "Ready Motor" (patented 1872) compressed air in one
cylinder, burned fuel in it at constant pressure, and expanded the gases in
another. The engine itself lost out to Otto's, but its cycle (also called
Joule's) is the cycle of every gas turbine and jet engine: isentropic
compression, heat addition at constant pressure, isentropic expansion and
heat rejection at constant pressure. With a calorically perfect gas the
efficiency depends only on the pressure ratio,
:math:`\eta = 1 - r_p^{(1-\gamma)/\gamma}`
(:func:`~engineeringkit.thermodynamics.brayton_efficiency`). The example
rebuilds that result from the four state temperatures.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import brayton_efficiency, isentropic_temperature_ratio, thermal_efficiency

# %%
# Four states from the first law
# ------------------------------
#
# Air at 300 K enters a compressor of pressure ratio 12; the turbine inlet is
# 1400 K. Heat in is :math:`c_p(T_3 - T_2)`, heat out :math:`c_p(T_4 - T_1)`.

cp, gamma, T1, T3, rp = 1005.0, 1.4, 300.0, 1400.0, 12.0
tau = isentropic_temperature_ratio(rp, gamma)
T2, T4 = T1 * tau, T3 / tau
q_in, q_out = cp * (T3 - T2), cp * (T4 - T1)
eta = thermal_efficiency(q_in, q_out)
print(f"T2 = {T2:.0f} K, T4 = {T4:.0f} K; net work {(q_in - q_out) / 1e3:.0f} kJ/kg")
print(f"efficiency from the states {eta:.4f}, from the formula {brayton_efficiency(rp, gamma):.4f}")
assert np.isclose(eta, brayton_efficiency(rp, gamma))
assert np.isclose(eta, 1 - T1 / T2)  # the efficiency is that of the compressor alone

# %%
# Efficiency against pressure ratio
# ---------------------------------
#
# A monatomic gas (:math:`\gamma = 5/3`), as in closed-cycle helium turbines,
# reaches a given efficiency at a lower pressure ratio than air. The
# efficiency keeps rising with pressure ratio, but the net work per kilogram
# does not: compressor work grows and, with the turbine inlet fixed, less heat
# can be added.

ratios = np.linspace(1, 40, 200)
fig, ax = plt.subplots(figsize=(6, 3.5))
for g, name in ((1.4, "air, γ = 1.4"), (5 / 3, "helium, γ = 5/3")):
    ax.plot(ratios, [brayton_efficiency(r, g) for r in ratios], label=name)
ax.axvline(rp, color="gray", ls=":")
ax.set(xlabel="pressure ratio", ylabel="ideal efficiency", ylim=(0, 1), title="Ideal Brayton cycle")
ax.legend()
fig.tight_layout()
assert brayton_efficiency(1.0) == 0 and brayton_efficiency(rp, 5 / 3) > brayton_efficiency(rp, 1.4)

plt.show()
