r"""
Stirling's regenerator (1816)
=============================

The Reverend Robert Stirling patented a hot-air engine in 1816 whose novelty
was the "economiser", a mass of metal through which the working air passed
on its way between the hot and cold spaces. Air cooling at constant volume
left its heat in the regenerator, and the returning air picked it up again.
Only the two isothermal heat exchanges then crossed the engine's boundary,
and the ideal cycle reaches Carnot's efficiency. Without the regenerator it
falls far short
(:func:`~engineeringkit.thermodynamics.stirling_efficiency`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import carnot_efficiency, isentropic_pressure_volume, stirling_efficiency
from engineeringkit.thermodynamics.visualizers import plot_pv_diagram

# %%
# The cycle on the p-v diagram
# ----------------------------
#
# One kilogram of air between 900 K and 300 K, compression ratio 4. The two
# isotherms are :math:`pv = RT`; the two constant-volume processes are the
# heat the regenerator stores and returns.

R, gamma, Th, Tc, r = 287.0, 1.4, 900.0, 300.0, 4.0
v_min = 0.1
v_max = r * v_min
v = np.linspace(v_min, v_max, 50)
ax = plot_pv_diagram(np.r_[v, v[::-1]], np.r_[R * Th / v, R * Tc / v[::-1]], label="Stirling")
ax.plot(v, isentropic_pressure_volume(R * Tc / v_max, v_max, v, gamma) / 1e3, ":", label="adiabat, for comparison")
ax.set(xlabel="specific volume [m³/kg]", title="Stirling cycle between 900 K and 300 K")
ax.legend()

w = R * (Th - Tc) * np.log(r)
q_regen = R / (gamma - 1) * (Th - Tc)
print(f"net work {w / 1e3:.1f} kJ/kg; heat stored and returned by the regenerator {q_regen / 1e3:.1f} kJ/kg")
for eps in (1.0, 0.9, 0.5, 0.0):
    print(f"regenerator effectiveness {eps:.1f}: efficiency {stirling_efficiency(Th, Tc, r, gamma, eps):.3f}")
print(f"Carnot {carnot_efficiency(Th, Tc):.3f}")
assert np.isclose(stirling_efficiency(Th, Tc, r), carnot_efficiency(Th, Tc))

# %%
# Why the regenerator matters
# ---------------------------
#
# The constant-volume heat is larger than the net work, so a poor regenerator
# throws most of the advantage away. A higher compression ratio helps only
# when the regenerator is imperfect.

eps = np.linspace(0, 1, 101)
fig, ax = plt.subplots(figsize=(6, 3.5))
for ratio in (2.0, 4.0, 10.0):
    ax.plot(eps, stirling_efficiency(Th, Tc, ratio, gamma, eps), label=f"r = {ratio:g}")
ax.axhline(carnot_efficiency(Th, Tc), color="gray", ls=":", label="Carnot")
ax.set(xlabel="regenerator effectiveness", ylabel="efficiency", title="Stirling engine, 900 K and 300 K")
ax.legend()
fig.tight_layout()
assert q_regen > w and stirling_efficiency(Th, Tc, r, gamma, 0.0) < 0.5 * carnot_efficiency(Th, Tc)

plt.show()
