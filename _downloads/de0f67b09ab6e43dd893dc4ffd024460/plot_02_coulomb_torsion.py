r"""
Coulomb's torsion
=================

Charles-Augustin Coulomb (1784) hung a weight on a thin wire, twisted it and
timed the oscillations. The restoring torque was proportional to the angle
of twist, to the fourth power of the wire's diameter and inversely to its
length, which in modern form is :math:`\phi = TL/(GJ)` with
:math:`J = \pi d^4/32`
(:func:`~engineeringkit.structural_engineering.torsion_angle`). The same
torsion balance later weighed the electric force and, in Cavendish's hands,
the Earth. The example simulates Coulomb's pendulum and recovers the shear
modulus from the period, then shows the shear stress across the wire
(:func:`~engineeringkit.structural_engineering.torsion_shear_stress`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from engineeringkit.structural_engineering import torsion_angle, torsion_shear_stress

# %%
# A torsion pendulum
# ------------------

L, d, G = 0.8, 0.5e-3, 40e9  # brass wire
J = np.pi * d**4 / 32
k = 1.0 / torsion_angle(1.0, L, G, J)  # torque per radian [N m/rad]
I_disk = 0.5 * 0.2 * 0.05**2  # 200 g disk of radius 50 mm [kg m^2]
sol = solve_ivp(lambda t, s: [s[1], -k * s[0] / I_disk], (0, 120), [0.5, 0.0], max_step=0.05, rtol=1e-9, atol=1e-12)
crossings = sol.t[1:][np.diff(np.sign(sol.y[0])) > 0]
period = float(np.mean(np.diff(crossings)))
G_measured = 4 * np.pi**2 * I_disk * L / (period**2 * J)
print(f"period {period:.2f} s, recovered G = {G_measured / 1e9:.2f} GPa")
assert np.isclose(G_measured, G, rtol=1e-3)

# %%
# Coulomb's fourth-power law
# --------------------------

diameters = np.array([0.25e-3, 0.5e-3, 1e-3])
torques = 0.1 / torsion_angle(1.0, L, G, np.pi * diameters**4 / 32)  # torque for 0.1 rad
print("torque ratios for doubling the diameter:", np.round(torques[1:] / torques[:-1], 3))
assert np.allclose(torques[1:] / torques[:-1], 16)

r = np.linspace(0, d / 2, 50)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
axes[0].plot(sol.t, sol.y[0])
axes[0].set(xlabel="time [s]", ylabel="twist [rad]", title="Torsion pendulum")
axes[1].plot(1e3 * r, torsion_shear_stress(k * 0.5, r, J) / 1e6)
axes[1].set(xlabel="radius [mm]", ylabel="shear stress [MPa]", title="Shear across the wire at 0.5 rad")
fig.tight_layout()

plt.show()
