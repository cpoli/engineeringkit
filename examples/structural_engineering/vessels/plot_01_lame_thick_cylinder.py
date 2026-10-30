r"""
Lamé's thick-walled cylinder (1833)
===================================

The hoop stress :math:`pr/t` of a boiler shell
(:func:`~engineeringkit.structural_engineering.thin_wall_vessel_stresses`)
assumes the stress is uniform through the wall. Gabriel Lamé and Émile
Clapeyron, in their memoir on the equilibrium of elastic solids (written
in Russia, published in Paris in 1833), solved the hollow cylinder exactly:
:math:`\sigma_r = A - B/r^2` and :math:`\sigma_\theta = A + B/r^2`
(:func:`~engineeringkit.structural_engineering.thick_wall_cylinder_stresses`).
The hoop stress is highest at the bore, and no thickness of wall can carry
an internal pressure above the material's strength. The example plots the
stresses through a gun-barrel wall and finds the radius ratio at which the
thin-wall formula is 5% low.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from engineeringkit.structural_engineering import thick_wall_cylinder_stresses, thin_wall_vessel_stresses

# %%
# Stresses through a thick wall
# -----------------------------

ri, ro, p = 0.05, 0.1, 200e6
r = np.linspace(ri, ro, 201)
radial, hoop, axial = thick_wall_cylinder_stresses(r, ri, ro, p)
print(f"hoop stress {hoop[0] / 1e6:.0f} MPa at the bore, {hoop[-1] / 1e6:.0f} MPa outside")
print(f"thin-wall formula {p * ri / (ro - ri) / 1e6:.0f} MPa")
assert np.isclose(hoop[0], p * (ro**2 + ri**2) / (ro**2 - ri**2))
assert np.isclose(np.trapezoid(hoop, r), p * ri)  # the wall carries the bursting force p * 2 ri per unit length

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(r * 1e3, hoop / 1e6, label="hoop")
ax.plot(r * 1e3, radial / 1e6, label="radial")
ax.plot(r * 1e3, axial / 1e6, label="axial (closed ends)")
ax.set(xlabel="radius [mm]", ylabel="stress [MPa]", title=f"Lamé cylinder, p = {p / 1e6:.0f} MPa")
ax.legend()
fig.tight_layout()

# %%
# When is a wall thin?
# --------------------


def error(ratio, t=1.0):
    """Shortfall of the thin-wall hoop stress against Lamé's peak, for inner radius / thickness."""
    exact = thick_wall_cylinder_stresses(ratio * t, ratio * t, (ratio + 1) * t, 1.0)[1]
    return 1 - thin_wall_vessel_stresses(1.0, ratio * t, t)[0] / exact


ratios = np.linspace(1, 40, 200)
limit = brentq(lambda k: error(k) - 0.05, 1, 40)
print(f"the thin-wall hoop stress is 5% low at r/t = {limit:.1f}")
assert 9 < limit < 11

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(ratios, 100 * np.array([error(k) for k in ratios]))
ax.axvline(limit, color="gray", ls="--")
ax.set(xlabel="inner radius / wall thickness", ylabel="thin-wall error [%]")
fig.tight_layout()

plt.show()
