"""
Poiseuille profile and volume flow
==================================

Integrate the laminar parabolic velocity profile over a circular cross section
and compare it with the volume-flow helper. The radius dependence is fourth power.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

from engineeringkit.fluid_engineering import poiseuille_flow

radius, length, mu, dp = 0.005, 1, 0.001, 10
r = np.linspace(0, radius, 10001)
v = dp * (radius**2 - r**2) / (4 * mu * length)
Q = poiseuille_flow(dp, radius, length, mu)
assert np.isclose(trapezoid(2 * np.pi * r * v, r), Q, rtol=1e-7)
assert np.isclose(poiseuille_flow(dp, 2 * radius, length, mu), 16 * Q)
fig, ax = plt.subplots()
ax.plot(r * 1000, v)
ax.set(xlabel="radial position [mm]", ylabel="axial speed [m/s]")

plt.show()
