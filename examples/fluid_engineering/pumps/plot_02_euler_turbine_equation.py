r"""
Euler's turbine equation (1754)
===============================

Leonhard Euler, studying Segner's reaction water wheel at the Berlin
academy, wrote the torque of a rotor as the change in the angular momentum
of the water passing through it. Multiplied by the speed, it gives the head
a pump adds or a turbine takes,
:math:`H = (u_2V_{t2} - u_1V_{t1})/g`
(:func:`~engineeringkit.fluid_engineering.euler_turbomachine_head`). This
example draws the ideal head of a centrifugal impeller for backward,
radial and forward blades, and shows that the head scales with the square
of the speed, the affinity law
(:func:`~engineeringkit.fluid_engineering.affinity_laws`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import affinity_laws, euler_turbomachine_head

# %%
# The outlet velocity triangle
# ----------------------------
#
# Water leaves an impeller of radius :math:`r_2` and width :math:`b_2`
# with radial speed :math:`V_{r2} = Q/(2\pi r_2 b_2)`. Relative to the
# blade it leaves at the blade angle :math:`\beta_2`, so its swirl is
# :math:`V_{t2} = u_2 - V_{r2}\cot\beta_2`. With no swirl at the inlet,
# the ideal head falls linearly with flow for backward-curved blades.

r2, b2, rpm = 0.15, 0.02, 1750
omega = rpm * 2 * np.pi / 60
u2 = omega * r2
Q = np.linspace(0, 0.08, 50)
V_r2 = Q / (2 * np.pi * r2 * b2)

fig, ax = plt.subplots(figsize=(6, 4))
heads = {}
for label, beta in [("backward, β₂ = 25°", 25), ("radial, β₂ = 90°", 90), ("forward, β₂ = 140°", 140)]:
    heads[beta] = euler_turbomachine_head(u2, u2 - V_r2 / np.tan(np.radians(beta)))
    ax.plot(Q * 1000, heads[beta], label=label)
ax.set(xlabel="flow [L/s]", ylabel="ideal head [m]", title=f"Euler head of an impeller at {rpm} rpm")
ax.legend()
fig.tight_layout()
print(f"shut-off head u2^2/g = {heads[25][0]:.1f} m for every blade angle")

# %%
# Backward blades are used in practice: their head falls with flow, so the
# pump is stable on its system curve and its power does not run away.
#
# The head at twice the speed
# ---------------------------
#
# At the same flow coefficient, the velocity triangle keeps its shape and
# every speed doubles, so the head goes up four times.

Q_bep, H_bep = 0.04, float(np.interp(0.04, Q, heads[25]))
Q2, H2, _ = affinity_laws(Q_bep, H_bep, 0.0, 2.0)
V_r2_fast = Q2 / (2 * np.pi * r2 * b2)
H2_euler = euler_turbomachine_head(2 * u2, 2 * u2 - V_r2_fast / np.tan(np.radians(25)))
print(f"at 2x speed: affinity law {H2:.1f} m, Euler {H2_euler:.1f} m")

assert np.allclose([heads[25][0], heads[90][0], heads[140][0]], u2**2 / 9.80665)
assert np.allclose(heads[90], heads[90][0]) and np.all(np.diff(heads[25]) < 0) and np.all(np.diff(heads[140]) > 0)
assert np.isclose(H2, H2_euler)

plt.show()
