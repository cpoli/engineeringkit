"""
Cooling and the lumped time constant
====================================

A small copper sphere has a small Biot number in the chosen convection field.
Its excess temperature should decrease to 1/e of its initial value at one
lumped time constant. All temperatures are absolute.
"""

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.heat_transfer import biot_number, lumped_temperature

radius, h, k, rho, cp = 0.005, 10, 400, 8960, 385
A, V = 4 * np.pi * radius**2, 4 * np.pi * radius**3 / 3
assert biot_number(h, V / A, k) < 0.1
tau = rho * cp * V / (h * A)
t = np.linspace(0, 5 * tau, 101)
T = lumped_temperature(373, 293, h, A, rho, cp, V, t)
assert np.isclose(lumped_temperature(373, 293, h, A, rho, cp, V, tau), 293 + 80 / np.e)
fig, ax = plt.subplots()
ax.plot(t / tau, (T - 293) / 80)
ax.set(xlabel="time / time constant", ylabel="normalized excess temperature")

plt.show()
