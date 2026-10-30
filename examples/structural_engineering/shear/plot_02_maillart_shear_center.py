r"""
Maillart's shear centre (1921)
==============================

Carl von Bach's tests on channel beams at Stuttgart showed them
twisting under loads applied through the web, and bending sideways under
loads in the plane of the web. Robert Maillart (1921) explained both. In a
section without two axes of symmetry the neutral axis is not perpendicular
to the plane of the loads unless the loads act along a principal axis
(:func:`~engineeringkit.structural_engineering.unsymmetrical_bending_stress`,
:func:`~engineeringkit.structural_engineering.neutral_axis_angle`); and the
shear flows in the flanges of a channel form a couple, so the load must pass
through a point outside the web, the *shear centre*
(:func:`~engineeringkit.structural_engineering.channel_shear_center`), to
bend the beam without twisting it. The example computes both for a channel
and an angle, built from rectangles by
:func:`~engineeringkit.structural_engineering.channel_section` and
:func:`~engineeringkit.structural_engineering.angle_section`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import (
    angle_section,
    channel_section,
    channel_shear_center,
    neutral_axis_angle,
    unsymmetrical_bending_stress,
)

# %%
# The channel: centroid and shear centre on opposite sides of the web
# --------------------------------------------------------------------

d, bf, tf, tw = 0.2, 0.075, 0.01, 0.006
channel = channel_section(d, bf, tf, tw)
h, b = d - tf, bf - tw / 2  # centreline dimensions
e = channel_shear_center(h, b, tf, tw)
print(f"centroid {channel.z_centroid * 1e3:.1f} mm in front of the web's back, shear centre {e * 1e3:.1f} mm behind its centreline")

# The flange shear flows rise linearly to V b h tf / (2 I) at the web; their couple F h must equal V e.
I_thin = tw * h**3 / 12 + 2 * b * tf * (h / 2) ** 2
F = b**2 * h * tf / (4 * I_thin)
assert np.isclose(F * h, e)

fig, ax = channel.plot()
ax.plot(tw / 2 - e, d / 2, "x", color="k", ms=10, label="shear centre")
ax.plot(channel.z_centroid, channel.y_centroid, "+", color="C3", ms=12, label="centroid")
ax.legend(loc="upper right")
ax.set_title("A channel must be loaded through its shear centre")

# %%
# The angle: a vertical load does not bend it vertically
# ------------------------------------------------------

angle = angle_section(0.15, 0.1, 0.012)
I1, I2, theta = angle.principal_moments()
Mz = -8e3  # hogging moment from a vertical load on a cantilever [N m]
alpha = neutral_axis_angle(0.0, Mz, angle.Iy, angle.Iz, angle.Iyz)
print(f"Iz {angle.Iz * 1e8:.0f} cm^4, Iy {angle.Iy * 1e8:.0f} cm^4, Iyz {angle.Iyz * 1e8:.0f} cm^4")
print(f"principal axes at {np.degrees(theta):.1f} deg, neutral axis at {np.degrees(alpha):.1f} deg (not 0)")
corners = np.array([[0, 0], [0, 0.15], [0.012, 0.15], [0.1, 0.012], [0.1, 0]]) - [angle.z_centroid, angle.y_centroid]
sigma = unsymmetrical_bending_stress(0.0, Mz, corners[:, 1], corners[:, 0], angle.Iy, angle.Iz, angle.Iyz)
for (z, y), stress in zip(corners, sigma):
    print(f"  corner ({z * 1e3:+.0f}, {y * 1e3:+.0f}) mm: {stress / 1e6:+.1f} MPa")
assert np.isclose(unsymmetrical_bending_stress(0.0, Mz, np.sin(alpha), np.cos(alpha), angle.Iy, angle.Iz, angle.Iyz), 0, atol=1e-3)
assert abs(alpha) > np.radians(5)

fig, ax = angle.plot()
t = np.linspace(-0.08, 0.08, 2)
ax.plot(angle.z_centroid + t * np.cos(alpha), angle.y_centroid + t * np.sin(alpha), "--", color="C3", label="neutral axis")
ax.plot(angle.z_centroid + t * np.cos(theta), angle.y_centroid + t * np.sin(theta), ":", color="k", label="major principal axis")
ax.legend(loc="upper right")
ax.set_title("Angle under a vertical load")

plt.show()
