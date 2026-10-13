"""
RLC current after a voltage step
================================

Compute an inspectable engineering model using SI inputs.
See the linked domain guide for assumptions and benchmark checks.
"""

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import SeriesRLC

t = np.linspace(0, 8, 801)
fig, ax = plt.subplots()
for resistance in [1, 2, 4]:
    result = SeriesRLC(resistance, 1, 0.25).step_response(t)
    ax.plot(result.t, result.response, label=f"R={resistance} ohm")
ax.set(xlabel="time [s]", ylabel="current [A] for a 1 V step")
ax.legend()

plt.show()
