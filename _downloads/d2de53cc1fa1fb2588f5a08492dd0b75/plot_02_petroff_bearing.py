r"""
Petroff's law of bearing friction (1883)
========================================

Railway axle bearings of the 1880s ran on oil, and their friction did not
behave like Coulomb's: it grew with speed and with the oil's viscosity, and
it fell as the load rose. Nikolai Petroff (1883) explained this by treating
the oil film as a thin layer sheared between the journal and its sleeve.
With the journal centred, the shear stress is :math:`\mu\omega r/c` all
round, so the friction torque is :math:`T = 2\pi\mu\omega r^3 L/c`
(:func:`~engineeringkit.mechanical_engineering.petroff_friction_torque`) and
the friction coefficient follows
:math:`f = 2\pi^2 (\mu N/P)(r/c)`. The example checks Petroff's thin-film
torque against the exact flow between concentric cylinders.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import petroff_friction_torque

# %%
# Thin film against exact Couette flow
# ------------------------------------
#
# Between cylinders of radii :math:`r` and :math:`r + c` the exact viscous
# torque is :math:`4\pi\mu L\omega r^2 (r + c)^2 / ((r + c)^2 - r^2)`. The
# relative difference from Petroff's formula is about :math:`3c/2r`, so the
# thin-film model is good to a fraction of a percent at the usual
# :math:`c/r \approx 10^{-3}`.

viscosity, r, L, omega = 0.03, 0.025, 0.05, 2 * np.pi * 25
ratios = np.geomspace(1e-4, 0.3, 40)
exact = np.array([4 * np.pi * viscosity * L * omega * r**2 * (r + c) ** 2 / ((r + c) ** 2 - r**2) for c in ratios * r])
petroff = petroff_friction_torque(viscosity, r, L, ratios * r, omega)
error = np.abs(petroff / exact - 1)
print(f"c/r = 1e-3: Petroff {petroff[np.searchsorted(ratios, 1e-3)]:.5f} N m, exact {exact[np.searchsorted(ratios, 1e-3)]:.5f} N m")
assert np.all(np.diff(error) > 0) and np.allclose(error[:10], 1.5 * ratios[:10], rtol=1e-2)

# %%
# The Petroff line
# ----------------
#
# Divide the torque by load times radius. Plotted against the bearing
# characteristic :math:`\mu N/P`, every load and speed falls on one line of
# slope :math:`2\pi^2 r/c`.

c = 25e-6
loads = np.array([500.0, 1000.0, 2000.0])
speeds = np.linspace(5, 60, 12)
fig, ax = plt.subplots(figsize=(6, 3.8))
for W in loads:
    P = W / (2 * r * L)
    f = petroff_friction_torque(viscosity, r, L, c, 2 * np.pi * speeds) / (W * r)
    characteristic = viscosity * speeds / P
    ax.plot(characteristic, f, "o", label=f"W = {W:.0f} N")
    assert np.allclose(f, 2 * np.pi**2 * characteristic * r / c)
line = np.linspace(0, ax.get_xlim()[1], 2)
ax.plot(line, 2 * np.pi**2 * line * r / c, "k--", lw=1, label="2 pi^2 (r/c) mu N / P")
ax.set(xlabel="bearing characteristic mu N / P [-]", ylabel="friction coefficient f", title="Petroff's law")
ax.legend()
fig.tight_layout()

plt.show()
