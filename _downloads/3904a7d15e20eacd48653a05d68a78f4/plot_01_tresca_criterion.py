r"""
Tresca's maximum-shear yield criterion
======================================

Henri Tresca (1864) punched and extruded lead, copper and iron through dies
at the Conservatoire des Arts et Métiers and found that metals flow like a
viscous fluid once the largest shear stress reaches a constant. With
principal stresses :math:`\sigma_1 \ge \sigma_2 \ge \sigma_3` the criterion
is :math:`\sigma_1 - \sigma_3 = \sigma_y`
(:func:`~engineeringkit.materials.tresca_equivalent`). It is the first
criterion for yielding under combined stresses, and in the plane-stress
principal plane it draws a hexagon. The example applies it to a shaft in
bending and torsion and to biaxial stress states.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import factor_of_safety, tresca_equivalent
from engineeringkit.materials.visualizers import plot_yield_surfaces

# %%
# A shaft in bending and torsion
# ------------------------------
#
# At the surface of a shaft the bending stress :math:`\sigma` and the
# torsional shear :math:`\tau` give principal stresses
# :math:`\sigma/2 \pm \sqrt{(\sigma/2)^2 + \tau^2}` and zero. Tresca's
# equivalent stress is their spread, :math:`\sqrt{\sigma^2 + 4\tau^2}`.

sy = 350e6
sigma, tau = 120e6, 90e6
center, radius = sigma / 2, np.hypot(sigma / 2, tau)
principal = [center + radius, center - radius, 0.0]
equivalent = tresca_equivalent(principal)
print(f"principal stresses {np.round(np.array(principal) / 1e6, 1)} MPa")
print(f"Tresca equivalent {equivalent / 1e6:.1f} MPa, safety factor {factor_of_safety(sy, equivalent):.2f}")
assert np.isclose(equivalent, np.sqrt(sigma**2 + 4 * tau**2))

# %%
# Where the zero stress matters
# -----------------------------
#
# When both in-plane principal stresses have the same sign, the largest
# shear acts out of the plane and the zero third stress sets the spread.
# Equal biaxial tension yields at :math:`\sigma_y`, not at infinity.

states = np.array([[1.0, 1.0, 0.0], [1.0, 0.5, 0.0], [1.0, -1.0, 0.0], [0.5, -0.5, 0.0]]) * sy
for state, eq in zip(states, tresca_equivalent(states)):
    print(f"sigma_1, sigma_2 = {state[0] / sy:+.1f}, {state[1] / sy:+.1f} sigma_y -> Tresca {eq / sy:.2f} sigma_y")
np.testing.assert_allclose(tresca_equivalent(states) / sy, [1, 1, 2, 1])

ax = plot_yield_surfaces(sy)
ax.plot(states[:, 0] / sy, states[:, 1] / sy, "ko")
ax.set_title("Tresca's hexagon")
ax.figure.tight_layout()

plt.show()
