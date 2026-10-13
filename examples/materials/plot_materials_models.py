"""
Voigt and Reuss composite bounds
================================

Compute an inspectable engineering model using SI inputs.
See the linked domain guide for assumptions and benchmark checks.
"""

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import rule_of_mixtures_E_longitudinal, rule_of_mixtures_E_transverse

fraction = np.linspace(0, 1, 101)
upper = [rule_of_mixtures_E_longitudinal(200e9, 3e9, f) / 1e9 for f in fraction]
lower = [rule_of_mixtures_E_transverse(200e9, 3e9, f) / 1e9 for f in fraction]
fig, ax = plt.subplots()
ax.plot(fraction, upper, label="longitudinal / Voigt")
ax.plot(fraction, lower, label="transverse / Reuss")
ax.set(xlabel="fibre volume fraction", ylabel="modulus [GPa]")
ax.legend()

plt.show()
