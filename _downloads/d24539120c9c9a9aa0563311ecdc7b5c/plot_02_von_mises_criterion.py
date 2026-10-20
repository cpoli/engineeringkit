r"""
Von Mises's yield criterion
===========================

Tresca's hexagon has corners, which make the equations of plasticity
awkward. Richard von Mises (1913) proposed a smooth surface through the
corners: yield when the second invariant of the deviatoric stress reaches
:math:`\sigma_y^2/3`, in plane stress
:math:`\sigma_x^2 - \sigma_x\sigma_y + \sigma_y^2 + 3\tau_{xy}^2 = \sigma_y^2`
(:func:`~engineeringkit.materials.von_mises_plane_stress`). Huber (1904)
had proposed the same rule, and Hencky (1924) explained it as the elastic
energy of distortion. Taylor and Quinney's tension-torsion tests on thin
tubes (1931) fell closer to von Mises than to Tresca.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import tresca_equivalent, von_mises_plane_stress
from engineeringkit.materials.visualizers import plot_yield_surfaces

# %%
# Distortion energy
# -----------------
#
# The strain energy per volume splits into a volume-change part and a
# distortion part, :math:`U_d = J_2/2G`. Setting :math:`U_d` equal to its
# value in a tensile test at yield gives the von Mises criterion.

G = 80e9
rng = np.random.default_rng(1)
sx, sy_, txy = rng.normal(0, 200e6, (3, 5))
s = np.array([[sx, txy, 0 * sx], [txy, sy_, 0 * sx], [0 * sx, 0 * sx, 0 * sx]]).transpose(2, 0, 1)
deviator = s - np.trace(s, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
J2 = 0.5 * np.einsum("kij,kij->k", deviator, deviator)
U_d = J2 / (2 * G)
equivalent = von_mises_plane_stress(sx, sy_, txy)
assert np.allclose(U_d, equivalent**2 / (6 * G))

# %%
# Tension and torsion of a thin tube
# ----------------------------------
#
# Under tension :math:`\sigma` and shear :math:`\tau` the two criteria are
# ellipses, :math:`\sigma^2 + 4\tau^2 = \sigma_y^2` (Tresca) and
# :math:`\sigma^2 + 3\tau^2 = \sigma_y^2` (von Mises). In pure shear they
# disagree most: von Mises predicts yield at :math:`\sigma_y/\sqrt3`,
# 15% above Tresca's :math:`\sigma_y/2`.

yield_strength = 300e6
tau_mises = yield_strength / np.sqrt(3)
print(f"pure-shear yield: von Mises {tau_mises / 1e6:.1f} MPa, Tresca {yield_strength / 2e6:.1f} MPa")
assert np.isclose(von_mises_plane_stress(0, 0, tau_mises), yield_strength)

angle = np.linspace(0, np.pi / 2, 91)
sig, tau = np.cos(angle), np.sin(angle)
mises_scale = 1 / von_mises_plane_stress(sig, 0, tau)
center, radius = sig / 2, np.hypot(sig / 2, tau)
tresca_scale = 1 / tresca_equivalent(np.stack([center + radius, center - radius, 0 * sig], axis=-1))
ratio = mises_scale / tresca_scale
print(f"von Mises allows between {ratio.min():.3f} and {ratio.max():.3f} times Tresca's load")
assert np.all(ratio >= 1 - 1e-12) and np.isclose(ratio.max(), 2 / np.sqrt(3))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4))
ax1.plot(sig * tresca_scale, tau * tresca_scale, label="Tresca")
ax1.plot(sig * mises_scale, tau * mises_scale, label="von Mises")
ax1.set(xlabel="sigma / sigma_y", ylabel="tau / sigma_y", title="Tension-torsion yield locus", aspect="equal")
ax1.legend()
plot_yield_surfaces(yield_strength, ax=ax2)
ax2.set_title("Plane-stress principal plane")
fig.tight_layout()

plt.show()
