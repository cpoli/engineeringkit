r"""
Torricelli's law of efflux (1644)
=================================

Evangelista Torricelli, Galileo's successor at Florence, found that water
spouts from a small hole in a vessel at the speed a body gains falling
from the free surface, :math:`v = \sqrt{2gh}`, and published it in the
*Opera geometrica* of 1644. The jet narrows just past a sharp-edged hole
(the vena contracta), so the flow is a discharge coefficient of about
0.61 times the ideal one
(:func:`~engineeringkit.fluid_engineering.orifice_flow`). This example
drains a tank through an orifice
(:func:`~engineeringkit.fluid_engineering.tank_drain_time`) and checks
the square-root law of the level against a step-by-step integration.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

from engineeringkit.fluid_engineering import G, orifice_flow, tank_drain_time

# %%
# Speed of the jet
# ----------------
#
# Torricelli's speed is the free-fall speed from the surface, so a jet
# pointed upward would just reach the surface again.

h = 2.0
v = np.sqrt(2 * G * h)
print(f"jet speed from {h} m of head: {v:.2f} m/s, rising {v**2 / (2 * G):.2f} m")

# %%
# Draining a tank
# ---------------
#
# A water butt 0.6 m across with a 20 mm sharp-edged hole at the bottom.
# The level obeys :math:`A_t\,dh/dt = -C_d A_o\sqrt{2gh}`, so
# :math:`\sqrt h` falls linearly in time.

A_tank, A_hole, h0, cd = np.pi * 0.3**2, np.pi * 0.01**2, 1.0, 0.61
T = tank_drain_time(A_tank, A_hole, h0, discharge_coefficient=cd)
t = np.linspace(0, T, 200)
h_exact = (np.sqrt(h0) - t / T * np.sqrt(h0)) ** 2

sol = integrate.solve_ivp(
    lambda t, h: [-orifice_flow(A_hole, max(h[0], 0.0), cd) / A_tank], (0, T), [h0], t_eval=t[:-1], rtol=1e-10, atol=1e-12
)
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(t / 60, h_exact, label="Torricelli, closed form")
ax.plot(sol.t[::20] / 60, sol.y[0][::20], "o", label="step-by-step integration")
ax.set(xlabel="time [min]", ylabel="water level [m]", title="Draining a tank through an orifice")
ax.legend()
fig.tight_layout()
print(f"drain time {T / 60:.1f} min; the top half goes in {tank_drain_time(A_tank, A_hole, h0, h0 / 2, cd) / 60:.1f} min")

# %%
# The vena contracta
# ------------------
#
# An ideal jet (:math:`C_d = 1`) would empty the tank 1/0.61 times faster.

ideal = tank_drain_time(A_tank, A_hole, h0, discharge_coefficient=1.0)
print(f"ideal drain time {ideal / 60:.1f} min")

assert np.allclose(sol.y[0], h_exact[:-1], atol=1e-6)
assert np.isclose(T / ideal, 1 / cd)
assert np.isclose(orifice_flow(A_hole, h, 1.0), A_hole * v)

plt.show()
