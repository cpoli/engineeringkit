r"""
Bredt's thin-walled torsion (1896)
==================================

Coulomb's :math:`\tau = Tr/J` holds only for circular shafts. Rudolf Bredt
(1896) showed that in any closed thin-walled tube the shear flow
:math:`q = \tau t` is the same all round the wall, :math:`q = T/(2A_m)`
with :math:`A_m` the area enclosed by the wall's centreline
(:func:`~engineeringkit.structural_engineering.thin_walled_torsion_shear_stress`),
and that the twist is :math:`\phi = (TL/4A_m^2G)\oint ds/t`
(:func:`~engineeringkit.structural_engineering.thin_walled_torsion_angle`).
Aircraft fuselages, box girders and car chassis rails are designed with
these two formulas. The example checks them against Coulomb for a circular
tube, then cuts the tube along its length: the open section, with torsion
constant :math:`st^3/3`, is hundreds of times more flexible.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import thin_walled_torsion_angle, thin_walled_torsion_shear_stress, torsion_angle

# %%
# A closed circular tube agrees with Coulomb
# ------------------------------------------

T, L, G = 1.5e3, 2.0, 80e9
r, t = 0.04, 0.002
Am, s = np.pi * r**2, 2 * np.pi * r
J_closed = 2 * np.pi * r**3 * t
phi_bredt = thin_walled_torsion_angle(T, L, G, Am, s / t)
print(f"Bredt {np.degrees(phi_bredt):.3f} deg, Coulomb {np.degrees(torsion_angle(T, L, G, J_closed)):.3f} deg")
print(f"shear stress {thin_walled_torsion_shear_stress(T, Am, t) / 1e6:.1f} MPa")
assert np.isclose(phi_bredt, torsion_angle(T, L, G, J_closed))

# %%
# Slitting the tube
# -----------------
#
# Cut along its length, the tube can only resist torque by shear across the
# thickness, like a flat strip: :math:`J = st^3/3`.

ratio = r / t * np.linspace(0.25, 1.5, 50)
J_open = s * t**3 / 3
stiffness_ratio = 4 * Am**2 / (s / t) / J_open
print(f"closed / slit stiffness at r/t = {r / t:.0f}: {stiffness_ratio:.0f}")
assert np.isclose(stiffness_ratio, 3 * (r / t) ** 2)

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.loglog(ratio, 3 * ratio**2)
ax.set(xlabel="radius / wall thickness", ylabel="closed / slit torsional stiffness", title="Why tubes are closed")
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()

plt.show()
