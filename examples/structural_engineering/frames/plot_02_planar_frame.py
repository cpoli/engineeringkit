r"""
Livesley's plane frames by direct stiffness (1956)
==================================================

Moment distribution assumes the joints do not move sideways. A portal
frame under wind does sway, and hand methods needed an extra correction
for every storey. In 1956 R. K. Livesley programmed the stiffness method
for rigid-jointed frames on the Manchester Mark I. Each member gets a 6x6
matrix with axial and bending terms, and the sway comes out with the
rotations (:class:`~engineeringkit.structural_engineering.PlanarFrame`).
The example sways a fixed-base portal frame and checks the stiffness
result against its two limits, an infinitely stiff beam and an infinitely
flexible one, then draws the deflected shape under wind and gravity.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import PlanarFrame

# %%
# A fixed-base portal frame
# -------------------------

h, L, E = 4.0, 6.0, 200e9
A, Ic = 1e-2, 8e-5  # columns
nodes = [[0, 0], [0, h], [L, h], [L, 0]]
members = ((0, 1), (1, 2), (2, 3))
fixed = [0, 1, 2, 9, 10, 11]
H = 20e3  # wind at the beam level [N]
loads = np.zeros(12)
loads[3] = H

# %%
# Sway against beam stiffness
# ---------------------------
#
# With an infinitely stiff beam, each column bends in double curvature and
# the frame sways :math:`Hh^3/(24EI_c)`. With no beam stiffness, each
# column is a cantilever: :math:`Hh^3/(6EI_c)`. The axial areas are made
# large so that only bending counts.

ratios = np.logspace(-4, 4, 33)  # I_beam / I_column
sway = []
for ratio in ratios:
    frame = PlanarFrame(nodes, members, E, 1e3 * A, [Ic, ratio * Ic, Ic])
    sway.append(frame.solve(loads, fixed).displacements[3])
sway = np.array(sway)
stiff, flexible = H * h**3 / (24 * E * Ic), H * h**3 / (6 * E * Ic)
print(f"sway: rigid-beam limit {1e3 * stiff:.2f} mm, stiffness method {1e3 * sway[-1]:.2f} mm")
print(f"sway: pinned-beam limit {1e3 * flexible:.2f} mm, stiffness method {1e3 * sway[0]:.2f} mm")
assert np.isclose(sway[-1], stiff, rtol=1e-3) and np.isclose(sway[0], flexible, rtol=1e-3)
assert np.all(np.diff(sway) < 0)

# %%
# Wind and gravity together
# -------------------------

frame = PlanarFrame(nodes, members, E, A, [Ic, 2 * Ic, Ic])
result = frame.solve(loads, fixed, member_loads=[0.0, 15e3, 0.0])
summary = result.summary()
print(f"sway {1e3 * result.displacements[3]:.2f} mm, largest end moment {summary['max_abs_end_moment'] / 1e3:.1f} kN m")
print(f"axial forces in the columns: {np.round(result.axial_forces()[[0, 2]] / 1e3, 1)} kN; cond(K_ff) = {summary['condition_number']:.1e}")
r = result.reactions
assert np.isclose(r[0] + r[9], -H)
assert np.isclose(r[1] + r[10], 15e3 * L)
# moments about the left base: base moments, vertical reaction, wind and gravity
assert np.isclose(r[2] + r[11] + r[10] * L - H * h - 15e3 * L * L / 2, 0)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].semilogx(ratios, 1e3 * sway, label="stiffness method")
axes[0].axhline(1e3 * stiff, color="gray", ls="--", label=r"$Hh^3/24EI_c$")
axes[0].axhline(1e3 * flexible, color="gray", ls=":", label=r"$Hh^3/6EI_c$")
axes[0].set(xlabel=r"$I_{beam}/I_{column}$", ylabel="sway [mm]", title="Portal frame sway")
axes[0].legend()
result.plot(scale=100, ax=axes[1])
axes[1].set_title("Deflected shape under wind and gravity (x100)")
fig.tight_layout()

plt.show()
