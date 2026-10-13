"""
Beam deflection under a uniform load
====================================

Compute an inspectable engineering model using SI inputs.
See the linked domain guide for assumptions and benchmark checks.
"""

import matplotlib.pyplot as plt

from engineeringkit.structural_engineering import SimplySupportedBeam

result = SimplySupportedBeam(5, 200e9, 8e-6).solve_uniform_load(2000)
print(result.summary())
fig, axes = result.plot()
for ax, label in zip(axes, ["shear [N]", "moment [N m]", "deflection [m]"]):
    ax.set_ylabel(label)
fig.tight_layout()

plt.show()
