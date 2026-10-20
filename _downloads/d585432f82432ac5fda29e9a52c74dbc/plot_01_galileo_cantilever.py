r"""
Galileo's cantilever
====================

In the *Two New Sciences* (1638) Galileo asked how strong a cantilever is,
and answered that its breaking moment grows with the width :math:`b` and the
*square* of the depth :math:`h`, so a plank on edge is stronger than one laid
flat by the ratio :math:`h/b`. He assumed the whole section pulls with the
same stress and pivots about its bottom edge, which gives
:math:`M = \sigma bh^2/2`. Linear elasticity (Navier, 1826) puts the pivot at
the centroid with a linear stress profile, :math:`M = \sigma bh^2/6`
(:func:`~engineeringkit.structural_engineering.rectangular_section_modulus`):
Galileo's scaling law was right, his constant three times too large.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson

from engineeringkit.structural_engineering import rectangular_second_moment, rectangular_section_modulus

# %%
# A plank on edge and laid flat
# -----------------------------

b, h, sigma = 0.05, 0.20, 10e6  # 50 x 200 mm timber, 10 MPa
for name, (width, depth) in {"on edge": (b, h), "flat": (h, b)}.items():
    Z = rectangular_section_modulus(width, depth)
    print(f"{name:8s}: Galileo M = {sigma * width * depth**2 / 2 / 1e3:5.1f} kN m, elastic M = {sigma * Z / 1e3:5.2f} kN m")
ratio = rectangular_section_modulus(b, h) / rectangular_section_modulus(h, b)
stiffness_ratio = rectangular_second_moment(b, h) / rectangular_second_moment(h, b)
print(f"strength ratio edge/flat = {ratio:.1f} (h/b = {h / b:.1f}), stiffness ratio = {stiffness_ratio:.0f}")
assert np.isclose(ratio, h / b) and np.isclose(stiffness_ratio, (h / b) ** 2)

# %%
# Where Galileo's factor of three comes from
# ------------------------------------------
#
# Integrate each assumed stress profile over the depth to get the moment.

y = np.linspace(-h / 2, h / 2, 2001)
galileo = np.full_like(y, sigma)  # uniform tension, pivot at the bottom edge
elastic = sigma * y / (h / 2)  # linear, pivot at the centroid
M_galileo = simpson(galileo * b * (y + h / 2), x=y)
M_elastic = simpson(elastic * b * y, x=y)
print(f"Galileo / elastic = {M_galileo / M_elastic:.3f}")
assert np.isclose(M_galileo / M_elastic, 3, rtol=1e-4)
assert np.isclose(M_elastic, sigma * rectangular_section_modulus(b, h), rtol=1e-4)

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(galileo / 1e6, y * 1e3, label="Galileo: uniform tension")
ax.plot(elastic / 1e6, y * 1e3, label="Navier: linear about centroid")
ax.axvline(0, color="gray", lw=0.5)
ax.set(xlabel="stress at failure [MPa]", ylabel="height in section [mm]", title="Stress across the root section")
ax.legend()
fig.tight_layout()

plt.show()
