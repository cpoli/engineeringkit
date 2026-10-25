r"""
Venturi's tube (1797)
=====================

Giovanni Battista Venturi reported in 1797 that water flowing through a
contraction drops in pressure at the narrowest point and recovers most of
it in a gentle expansion. Clemens Herschel turned the effect into the
Venturi water meter in 1887: continuity and Bernoulli's equation give the
flow from the pressure drop between inlet and throat,
:math:`Q = C_d A_2\sqrt{2\Delta p/[\rho(1 - (A_2/A_1)^2)]}`
(:func:`~engineeringkit.fluid_engineering.venturi_flow`). This example
draws the pressure along a meter and reads the flow back from it.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import bernoulli_pressure, continuity_velocity, venturi_flow

# %%
# Pressure along the meter
# ------------------------
#
# A 150 mm pipe narrows to a 75 mm throat over 0.2 m and widens again over
# 0.6 m, carrying 30 L/s of water. Without losses the pressure follows the
# area.

rho, Q, p_in = 998.0, 0.030, 200e3
x = np.linspace(0, 1.0, 201)
d = np.interp(x, [0, 0.1, 0.3, 0.35, 0.95, 1.0], [0.15, 0.15, 0.075, 0.075, 0.15, 0.15])
A = np.pi * d**2 / 4
v = continuity_velocity(A[0], Q / A[0], A)
p = bernoulli_pressure(p_in, rho, v[0], v)

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(x, p / 1e3)
ax.set(xlabel="position [m]", ylabel="pressure [kPa]", title="Ideal pressure along a Venturi meter")
fig.tight_layout()

dp = p[0] - p.min()
print(f"throat speed {v.max():.2f} m/s, pressure drop {dp / 1e3:.1f} kPa")

# %%
# Reading the flow
# ----------------
#
# The meter inverts the pressure drop. A real tube reads about 2% low
# because of friction, which the discharge coefficient corrects.

A1, A2 = A[0], A.min()
print(f"ideal meter: {1e3 * venturi_flow(dp, rho, A1, A2):.2f} L/s")
print(f"with C_d = 0.98: {1e3 * venturi_flow(dp, rho, A1, A2, 0.98):.2f} L/s")

assert np.isclose(venturi_flow(dp, rho, A1, A2), Q)
assert np.isclose(p[-1], p_in) and np.isclose(v.max(), 4 * v[0])

plt.show()
