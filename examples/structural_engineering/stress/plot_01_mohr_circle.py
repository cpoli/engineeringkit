r"""
Mohr's circle (1882)
====================

Otto Mohr (1882) showed that the normal and shear stresses on every plane
through a point, plotted against each other, fall on a single circle. Its
centre is the mean normal stress, its radius the largest in-plane shear
stress, and its ends on the axis are the principal stresses
(:func:`~engineeringkit.structural_engineering.mohr_circle`). Turning the
element by :math:`\theta` moves the point by :math:`2\theta` around the
circle (:func:`~engineeringkit.structural_engineering.stress_transformation`).
The example traces the circle by rotating an element, checks that the
principal stresses are the eigenvalues of the stress tensor, and finds the
plane of maximum shear.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import mohr_circle, stress_transformation
from engineeringkit.structural_engineering.visualizers import plot_mohr_circle

# %%
# A shaft under bending and torsion
# ---------------------------------
#
# At the surface of a shaft, bending gives :math:`\sigma_x` and torsion
# gives :math:`\tau_{xy}`; there is no stress across the shaft.

sx, sy, txy = 120e6, 0.0, 50e6
center, radius, s1, s2 = mohr_circle(sx, sy, txy)
print(f"centre {center / 1e6:.1f} MPa, radius {radius / 1e6:.1f} MPa, principal {s1 / 1e6:.1f} and {s2 / 1e6:.1f} MPa")
eigen = np.linalg.eigvalsh([[sx, txy], [txy, sy]])
assert np.allclose([s2, s1], eigen)

theta = np.linspace(0, np.pi, 181)
sxp, syp, txyp = stress_transformation(sx, sy, txy, theta)
assert np.allclose(np.hypot(sxp - center, txyp), radius)
assert np.allclose(sxp + syp, sx + sy)  # the first invariant

# %%
# Principal planes and the plane of maximum shear
# -----------------------------------------------

theta_p = 0.5 * np.arctan2(2 * txy, sx - sy)
on_principal = stress_transformation(sx, sy, txy, theta_p)
on_shear = stress_transformation(sx, sy, txy, theta_p + np.pi / 4)
print(f"principal plane at {np.degrees(theta_p):.1f} deg: shear {on_principal[2]:.1e} Pa")
print(f"max-shear plane at {np.degrees(theta_p) + 45:.1f} deg: shear {on_shear[2] / 1e6:.1f} MPa")
assert np.isclose(on_principal[0], s1) and abs(on_principal[2]) < 1e-6 * radius
assert np.isclose(abs(on_shear[2]), radius) and np.isclose(on_shear[0], center)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
plot_mohr_circle(sx, sy, txy, ax=axes[0])
axes[0].plot(sxp[::15], txyp[::15], ".", color="k", label="rotated element")
axes[0].set_title("Mohr's circle")
axes[1].plot(np.degrees(theta), sxp / 1e6, label=r"$\sigma_{x'}$")
axes[1].plot(np.degrees(theta), txyp / 1e6, label=r"$\tau_{x'y'}$")
axes[1].axvline(np.degrees(theta_p), color="gray", ls="--")
axes[1].set(xlabel="element rotation [deg]", ylabel="stress [MPa]", title="Stress transformation")
axes[1].legend()
fig.tight_layout()

plt.show()
