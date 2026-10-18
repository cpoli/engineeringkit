r"""
Pitot's tube
============

In 1732 Henri Pitot lowered a glass tube bent at a right angle into the
Seine, its mouth facing upstream, and watched the water inside rise above
the river surface. The rise is the velocity head :math:`v^2/(2g)`, and it
showed, against the belief of the day, that the current is faster near the
surface than at depth. Turned around, the stagnation rise gives the speed,
:math:`v = \sqrt{2\Delta p/\rho}`
(:func:`~engineeringkit.fluid_engineering.pitot_velocity`). This example
reads river and air speeds from a Pitot tube and finds the Mach number at
which the incompressible formula starts to mislead an airspeed indicator.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import G, dynamic_pressure, pitot_velocity, velocity_head

# %%
# Water rising in the tube
# ------------------------

for v in [0.5, 1.0, 2.0]:
    print(f"river at {v:.1f} m/s: water rises {100 * velocity_head(v):.1f} cm in the tube")

rise = 0.051
print(f"a {100 * rise:.1f} cm rise means {pitot_velocity(998.0 * G * rise, 998.0):.2f} m/s")

# %%
# Airspeed and compressibility
# ----------------------------
#
# In a gas the stagnation rise is
# :math:`p[(1 + \frac{\gamma - 1}{2}M^2)^{\gamma/(\gamma - 1)} - 1]`, larger
# than :math:`\tfrac12\rho v^2`. Feeding it into the incompressible formula
# overestimates the speed.

gamma, p, T, R = 1.4, 101325.0, 288.15, 287.05
rho, a = p / (R * T), np.sqrt(gamma * R * T)
M = np.linspace(0.01, 0.8, 80)
dp = p * ((1 + (gamma - 1) / 2 * M**2) ** (gamma / (gamma - 1)) - 1)
error = pitot_velocity(dp, rho) / (M * a) - 1

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(M, 100 * error)
ax.axhline(1, color="k", ls=":")
ax.set(xlabel="Mach number", ylabel="overestimate of speed [%]", title="Incompressible Pitot formula in air")
fig.tight_layout()

M_1pct = np.interp(0.01, error, M)
print(f"the error reaches 1% at Mach {M_1pct:.2f}")

# %%
# The tube inverts the dynamic pressure exactly, and stays within 1% in air
# below about Mach 0.28.

v = np.array([0.3, 3.0, 30.0])
assert np.allclose(pitot_velocity(dynamic_pressure(998.0, v), 998.0), v)
assert 0.25 < M_1pct < 0.3

plt.show()
