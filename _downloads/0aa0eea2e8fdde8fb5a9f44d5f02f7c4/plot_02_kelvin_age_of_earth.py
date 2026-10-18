r"""
Kelvin's age of the Earth (1862)
================================

In 1862 William Thomson, later Lord Kelvin, treated the Earth as a
semi-infinite solid that began at a uniform molten temperature and has
cooled from its surface ever since. The temperature below the surface is
then :math:`T_s + (T_i - T_s)\operatorname{erf}(x/2\sqrt{\alpha t})`
(:func:`~engineeringkit.heat_transfer.semi_infinite_temperature`), and the
surface gradient :math:`(T_i - T_s)/\sqrt{\pi\alpha t}` falls with time. From
the gradient measured in mines, about 1 °F per 50 ft, he found an age of
about 98 million years. This example redoes his arithmetic in SI units and
checks the gradient against the profile itself.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.heat_transfer import semi_infinite_temperature

# %%
# Kelvin's numbers
# ----------------
#
# An initial excess of 7000 °F, a gradient of 1 °F per 50 ft and a
# diffusivity of 400 ft² per year, converted to SI.

year = 365.25 * 86400
dT = 7000 * 5 / 9  # K
gradient = (5 / 9) / (50 * 0.3048)  # K/m
alpha = 400 * 0.3048**2 / year  # m^2/s
age = dT**2 / (np.pi * alpha * gradient**2)
print(f"diffusivity {alpha:.3e} m²/s, gradient {gradient * 1e3:.1f} mK/m")
print(f"age of the Earth: {age / year / 1e6:.0f} million years")
assert 95e6 < age / year < 100e6

# %%
# The gradient from the profile
# -----------------------------
#
# Differentiate the error-function profile numerically at the surface and
# compare with the closed form.

Ts, Ti = 273.0, 273.0 + dT
dx = 1.0
numeric = (semi_infinite_temperature(Ti, Ts, dx, alpha, age) - semi_infinite_temperature(Ti, Ts, 0.0, alpha, age)) / dx
print(f"surface gradient from the profile {numeric * 1e3:.3f} mK/m")
assert np.isclose(numeric, gradient, rtol=1e-4)

depth = np.linspace(0, 400e3, 400)
fig, ax = plt.subplots(figsize=(6, 4))
for myr in [1, 10, 98, 1000]:
    ax.plot(semi_infinite_temperature(Ti, Ts, depth, alpha, myr * 1e6 * year) - 273.0, depth / 1e3, label=f"{myr} Myr")
ax.invert_yaxis()
ax.set(xlabel="temperature [°C]", ylabel="depth [km]", title="Kelvin's cooling Earth")
ax.legend()
fig.tight_layout()

# %%
# Why he was wrong
# ----------------
#
# Radioactive heating (discovered in 1896) and, above all, convection in the
# mantle, argued by John Perry in 1895, keep the interior hot. The radiometric
# age, 4.54 billion years, would need a surface gradient about
# :math:`\sqrt{4540/98}` times smaller than the one measured.

print(f"the gradient for 4.54 Gyr would be {gradient * np.sqrt(age / (4.54e9 * year)) * 1e3:.1f} mK/m")

plt.show()
