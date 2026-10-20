r"""
Stevin's hydrostatic paradox
============================

In *De Beghinselen des Waterwichts* (1586) Simon Stevin showed that the
force of water on the bottom of a vessel depends only on the depth and the
area of the bottom, not on how much water the vessel holds:
:math:`p = \rho g h` (:func:`~engineeringkit.fluid_engineering.hydrostatic_pressure`).
He went on to find the force on a sloping or vertical wall, the first
problem of engineering hydrostatics. This example fills three vessels of
different shape to the same depth, then computes the thrust on a dam face
and where it acts (:func:`~engineeringkit.fluid_engineering.hydrostatic_force_plane`,
:func:`~engineeringkit.fluid_engineering.center_of_pressure_depth`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import (
    G,
    center_of_pressure_depth,
    hydrostatic_force_plane,
    hydrostatic_pressure,
    pressure_head,
)

# %%
# Three vessels, one base force
# -----------------------------
#
# Each vessel has a 0.01 m² base and holds water 0.5 m deep. A flared vessel
# holds five times the water of a narrow-necked one, yet the force on the
# base is the same; the side walls carry the difference.

rho, depth, base = 998.0, 0.5, 0.01
volumes = {"narrow neck": 0.002, "straight": 0.005, "flared": 0.010}
base_force = hydrostatic_pressure(rho, depth) * base
for name, V in volumes.items():
    print(f"{name:>11}: water weighs {rho * G * V:5.1f} N, force on base {base_force:.1f} N")

# %%
# A dam face
# ----------
#
# Water 30 m deep stands against a vertical dam 1 m wide. The pressure grows
# linearly with depth, so the thrust is that at mid-depth times the area and
# acts two thirds of the way down. A 30 m column is about 3 atmospheres.

H, width = 30.0, 1.0
F = hydrostatic_force_plane(rho, H / 2, H * width)
y_p = center_of_pressure_depth(H / 2, H * width, width * H**3 / 12)
print(f"thrust {F / 1e6:.2f} MN per metre, acting {y_p:.1f} m below the surface")
p_base = hydrostatic_pressure(rho, H)
print(f"pressure at the base {p_base / 1e3:.0f} kPa = {pressure_head(p_base, 13546.0):.2f} m of mercury")

y = np.linspace(0, H, 50)
fig, ax = plt.subplots(figsize=(5, 4))
ax.fill_betweenx(y, 0, hydrostatic_pressure(rho, y) / 1e3, alpha=0.4)
ax.axhline(y_p, color="k", ls="--", label=f"centre of pressure, {y_p:.0f} m")
ax.invert_yaxis()
ax.set(xlabel="pressure [kPa]", ylabel="depth [m]", title="Pressure on a dam face")
ax.legend()
fig.tight_layout()

# %%
# Stevin's result: the thrust on the vertical face is half the force the
# same area would carry lying flat at the full depth.

assert np.isclose(F, 0.5 * hydrostatic_pressure(rho, H) * H * width)
assert np.isclose(y_p, 2 * H / 3)
assert np.isclose(base_force, rho * G * depth * base)

plt.show()
