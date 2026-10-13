"""
Carnot engine and heat-pump limits
==================================

A reversible heat pump delivers absorbed heat plus supplied work.
Consequently its COP exceeds refrigerator COP by one. The temperature
difference must remain positive and absolute temperatures are used.
"""

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import carnot_cop_heat_pump, carnot_cop_refrigerator, carnot_efficiency

Tc = 280
Th = np.linspace(281, 450, 201)
engine = np.array([carnot_efficiency(t, Tc) for t in Th])
pump = np.array([carnot_cop_heat_pump(t, Tc) for t in Th])
fridge = np.array([carnot_cop_refrigerator(t, Tc) for t in Th])
np.testing.assert_allclose(pump - fridge, 1)
assert np.isclose(carnot_efficiency(560, Tc), 0.5)
fig, axes = plt.subplots(1, 2, figsize=(9, 3))
axes[0].plot(Th, engine)
axes[0].set(xlabel="hot reservoir [K]", ylabel="engine efficiency")
axes[1].plot(Th, pump, label="heat pump")
axes[1].plot(Th, fridge, label="refrigerator")
axes[1].set(xlabel="hot reservoir [K]", ylabel="reversible COP", ylim=(0, 30))
axes[1].legend()
fig.tight_layout()

plt.show()
