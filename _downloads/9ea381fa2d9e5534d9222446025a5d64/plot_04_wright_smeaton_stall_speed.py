r"""
The Wright brothers, Smeaton's coefficient and the stall speed
==============================================================

The nineteenth-century lift equation was :math:`L = k V^2 S C_L` with
Smeaton's coefficient :math:`k = 0.005` lb/(ft² mph²). After their 1901
glider lifted far less than predicted, Wilbur and Orville Wright measured
:math:`k \approx 0.0033` in their own wind tunnel; the modern value,
:math:`\tfrac12\rho` at sea level, is 0.00256. Since the minimum flying speed
:math:`V_s = \sqrt{2W/(\rho S C_{L,\max})}`
(:func:`~engineeringkit.aerodynamics.stall_speed`) scales as
:math:`1/\sqrt{k}`, Smeaton's value under-predicts the take-off speed of the
1903 Flyer by about 30%.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import stall_speed
from engineeringkit.units import convert

# %%
# Smeaton's coefficient as an implied air density
# -----------------------------------------------
#
# :math:`k = \tfrac12\rho`, so each historical value implies a density.

lbf_per_ft2 = convert(1.0, "psi", "Pa") / 144
mph = 0.44704
for label, k in [("Smeaton (1759)", 0.005), ("Wright tunnel (1901)", 0.0033), ("modern", 0.00256)]:
    rho_implied = 2 * k * lbf_per_ft2 / mph**2
    print(f"{label:22s} k = {k:.5f} -> rho = {rho_implied:.3f} kg/m^3")

# %%
# The 1903 Flyer's minimum speed
# ------------------------------
#
# About 340 kg flying mass, 47.4 m² of wing, and :math:`C_{L,\max}` taken
# between 0.6 and 1.4.

weight = 340 * 9.80665
S = 47.4
clmax = np.linspace(0.6, 1.4, 50)
fig, ax = plt.subplots(figsize=(6, 4))
for label, k in [("Smeaton k = 0.005", 0.005), ("Wright k = 0.0033", 0.0033), ("modern k = 0.00256", 0.00256)]:
    rho_implied = 2 * k * lbf_per_ft2 / mph**2
    ax.plot(clmax, stall_speed(weight, rho_implied, S, clmax), label=label)
ratio = float(stall_speed(weight, 1.225, S, 1.0) / stall_speed(weight, 2 * 0.005 * lbf_per_ft2 / mph**2, S, 1.0))
print(f"true / Smeaton-predicted stall speed: {ratio:.2f}")
ax.set(xlabel="$C_{L,\\max}$", ylabel="stall speed [m/s]", title="1903 Flyer minimum flying speed")
ax.legend()
fig.tight_layout()

plt.show()
