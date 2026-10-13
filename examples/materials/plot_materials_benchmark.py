"""
Crack size and fracture threshold
=================================

With a fixed tensile stress and geometry factor, stress intensity rises with
the square root of crack length. Inverting the toughness threshold recovers
the critical crack length. Small-scale yielding is assumed.
"""

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import critical_crack_length, stress_intensity

sigma, Y, Kic = 100e6, 1.12, 50e6
critical = critical_crack_length(Kic, sigma, Y)
a = np.linspace(1e-6, 1.5 * critical, 201)
K = np.array([stress_intensity(Y, sigma, ai) for ai in a])
assert np.isclose(stress_intensity(Y, sigma, critical), Kic)
fig, ax = plt.subplots()
ax.plot(a * 1000, K / 1e6)
ax.axhline(Kic / 1e6, color="black", linestyle="--", label="toughness")
ax.axvline(critical * 1000, color="gray", linestyle=":", label="critical crack")
ax.set(xlabel="crack length [mm]", ylabel="stress intensity [MPa sqrt(m)]")
ax.legend()

plt.show()
