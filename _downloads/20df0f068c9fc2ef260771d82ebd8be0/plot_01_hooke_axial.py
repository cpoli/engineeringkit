r"""
Hooke's law: ut tensio, sic vis
===============================

In 1676 Robert Hooke published an anagram, *ceiiinosssttuu*, and in 1678 its
solution, *ut tensio, sic vis*: "as the extension, so the force". Springs,
wires and wooden rods all stretched in proportion to the weight hung on
them. Written for a bar of length :math:`L`, area :math:`A` and modulus
:math:`E` (Young's, 1807), the extension is :math:`\delta = PL/(AE)`
(:func:`~engineeringkit.structural_engineering.axial_deformation`) and the
stored energy is the triangle under the load-extension line,
:math:`U = P\delta/2 = P^2L/(2AE)`
(:func:`~engineeringkit.structural_engineering.strain_energy_axial`). The
example repeats Hooke's experiment on a steel wire and checks both.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson

from engineeringkit.structural_engineering import axial_deformation, strain_energy_axial

# %%
# Weights on a wire
# -----------------

L, d, E = 2.0, 1e-3, 200e9  # a 2 m steel wire, 1 mm diameter
A = np.pi * d**2 / 4
loads = 9.81 * np.arange(0, 11) * 2.0  # 0 to 20 kg in 2 kg steps [N]
stretch = axial_deformation(loads, L, A, E)
for P, dl in list(zip(loads, stretch))[::5]:
    print(f"load {P:6.1f} N: extension {1e3 * dl:.3f} mm")
stiffness = np.polyfit(stretch, loads, 1)[0]
print(f"measured stiffness {stiffness / 1e3:.1f} kN/m = AE/L = {A * E / L / 1e3:.1f} kN/m")
assert np.isclose(stiffness, A * E / L)

# %%
# Energy is the area under the line
# ---------------------------------

P_max = loads[-1]
area = simpson(loads, x=stretch)
U = strain_energy_axial(P_max, L, A, E)
print(f"area under the load-extension line {1e3 * area:.3f} mJ, P^2 L/(2AE) = {1e3 * U:.3f} mJ")
assert np.isclose(area, U) and np.isclose(U, P_max * stretch[-1] / 2)

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(1e3 * stretch, loads, "o-")
ax.fill_between(1e3 * stretch, loads, alpha=0.2, label=f"U = {1e3 * U:.1f} mJ")
ax.set(xlabel="extension [mm]", ylabel="load [N]", title="Ut tensio, sic vis")
ax.legend()
fig.tight_layout()

plt.show()
