"""
Parallel and counterflow heat exchangers
========================================

Compute an inspectable engineering model using SI inputs.
See the linked domain guide for assumptions and benchmark checks.
"""

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.heat_transfer import exchanger_effectiveness_counterflow, exchanger_effectiveness_parallel

ntu = np.linspace(0, 5, 101)
fig, ax = plt.subplots()
for ratio in [0.5, 1]:
    counter = [exchanger_effectiveness_counterflow(n, ratio) for n in ntu]
    parallel = [exchanger_effectiveness_parallel(n, ratio) for n in ntu]
    ax.plot(ntu, counter, label=f"counterflow Cr={ratio}")
    ax.plot(ntu, parallel, "--", label=f"parallel Cr={ratio}")
ax.set(xlabel="NTU", ylabel="effectiveness", ylim=(0, 1))
ax.legend()

plt.show()
