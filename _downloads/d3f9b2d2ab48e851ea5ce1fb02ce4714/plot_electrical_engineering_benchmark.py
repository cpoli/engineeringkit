"""
RC charge and energy accounting
===============================

Charging an initially empty ideal capacitor from a constant voltage source
dissipates half the eventual source energy in the resistor. Numerical integration
checks energy conservation over a finite time interval.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

from engineeringkit.electrical_engineering import rc_charging_voltage, rc_time_constant

R, C, Vs = 1000, 1e-6, 5
tau = rc_time_constant(R, C)
t = np.linspace(0, 10 * tau, 10001)
v = rc_charging_voltage(Vs, R, C, t)
i = (Vs - v) / R
source = trapezoid(Vs * i, t)
loss = trapezoid(R * i**2, t)
stored = 0.5 * C * v[-1] ** 2
assert np.isclose(source, loss + stored, rtol=1e-6)
assert np.isclose(rc_charging_voltage(Vs, R, C, tau), Vs * (1 - np.exp(-1)))
fig, ax = plt.subplots()
ax.plot(t / tau, v, label="capacitor voltage")
ax.set(xlabel="time / RC", ylabel="voltage [V]")
ax.legend()

plt.show()
