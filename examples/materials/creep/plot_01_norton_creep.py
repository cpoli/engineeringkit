r"""
Norton's power-law creep (1929)
===============================

Steam turbines of the 1920s ran hot enough for their steel to stretch slowly
under constant load. Frederick Norton's *The Creep of Steel at High
Temperatures* (1929) found that the steady creep rate rises as a power of
stress, :math:`\dot\varepsilon = A\sigma^n`, with :math:`n` around 5. Its
temperature dependence is thermally activated,
:math:`e^{-Q/RT}`, with an activation energy close to that of
self-diffusion (Dorn, 1950s)
(:func:`~engineeringkit.materials.arrhenius_creep_rate`). The example
recovers :math:`n` and :math:`Q` from simulated creep tests.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import arrhenius_creep_rate
from engineeringkit.materials.utils import power_law_fit

# %%
# Creep tests at three temperatures
# ---------------------------------
#
# Minimum creep rates of a ferritic steel with :math:`n = 5` and
# :math:`Q = 280` kJ/mol, measured with 15% scatter at five stresses and
# three temperatures.

A, Q, R, n = 1.8e-32, 280e3, 8.314, 5.0
stresses = np.geomspace(60e6, 200e6, 5)
temperatures = np.array([850.0, 900.0, 950.0])
rng = np.random.default_rng(11)
rates = arrhenius_creep_rate(A, Q, R, temperatures[:, None], stresses[None, :], n) * rng.lognormal(0, 0.15, (3, 5))

# %%
# The stress exponent and the activation energy
# ---------------------------------------------
#
# At each temperature the slope of :math:`\log\dot\varepsilon` against
# :math:`\log\sigma` is :math:`n`. At fixed stress the slope of
# :math:`\ln\dot\varepsilon` against :math:`1/T` is :math:`-Q/R`.

exponents = [power_law_fit(stresses, row)[1] for row in rates]
slope = np.polyfit(1 / temperatures, np.log(rates).mean(axis=1), 1)[0]
print(f"stress exponents {np.round(exponents, 2).tolist()}; activation energy {-slope * R / 1e3:.0f} kJ/mol")
assert np.allclose(exponents, n, atol=0.4) and abs(-slope * R / Q - 1) < 0.05

life = 0.01 / arrhenius_creep_rate(A, Q, R, 850.0, 100e6, n)
print(f"time to 1% creep strain at 850 K and 100 MPa: {life / 3600:.0f} h")

fig, ax = plt.subplots(figsize=(6, 3.8))
for T, row in zip(temperatures, rates):
    ax.loglog(stresses / 1e6, row, "o", label=f"{T:.0f} K")
    ax.loglog(stresses / 1e6, arrhenius_creep_rate(A, Q, R, T, stresses, n), "k-", lw=0.8)
ax.set(xlabel="stress [MPa]", ylabel="minimum creep rate [1/s]", title="Norton's law")
ax.legend()
fig.tight_layout()

plt.show()
