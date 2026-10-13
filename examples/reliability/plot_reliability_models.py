"""
Weibull survival and shape
==========================

Compute an inspectable engineering model using SI inputs.
See the linked domain guide for assumptions and benchmark checks.
"""

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.reliability import weibull_reliability

t = np.linspace(0, 8, 201)
fig, ax = plt.subplots()
for beta in [0.5, 1, 2]:
    ax.plot(t, weibull_reliability(t, 3, beta), label=f"shape beta={beta}")
ax.set(xlabel="time [s]", ylabel="survival probability", ylim=(0, 1))
ax.legend()

plt.show()
