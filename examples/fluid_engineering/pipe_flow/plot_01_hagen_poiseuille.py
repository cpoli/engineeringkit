r"""
Hagen and Poiseuille: laminar pipe flow (1839-1840)
===================================================

Gotthilf Hagen in 1839, timing water through brass tubes, and Jean Léonard
Marie Poiseuille in 1840, studying blood flow in glass capillaries, found
independently that slow flow through a tube grows with the pressure drop
and with the fourth power of the radius. The law
:math:`Q = \pi R^4\Delta p/(8\mu L)`
(:func:`~engineeringkit.fluid_engineering.poiseuille_flow`) follows from a
parabolic velocity profile
(:func:`~engineeringkit.fluid_engineering.poiseuille_velocity`), and in
Darcy's form it is the laminar friction factor :math:`f = 64/Re`
(:func:`~engineeringkit.fluid_engineering.laminar_friction_factor`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

from engineeringkit.fluid_engineering import (
    G,
    darcy_head_loss,
    laminar_friction_factor,
    poiseuille_flow,
    poiseuille_velocity,
    reynolds_number,
)

# %%
# The parabolic profile
# ---------------------
#
# Water (:math:`\mu = 1.0\times10^{-3}` Pa s) in a 2 mm capillary, 0.5 m
# long, with 20 Pa across it.

rho, mu, R, L, dp = 998.0, 1.0e-3, 1.0e-3, 0.5, 20.0
r = np.linspace(-R, R, 101)
u = poiseuille_velocity(np.abs(r), dp, R, L, mu)
Q = poiseuille_flow(dp, R, L, mu)
v_mean = Q / (np.pi * R**2)
print(f"Q = {Q * 1e9:.2f} mm³/s, mean speed {1e3 * v_mean:.2f} mm/s, centreline {1e3 * u.max():.2f} mm/s")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.8))
ax1.plot(1e3 * u, 1e3 * r)
ax1.axvline(1e3 * v_mean, color="k", ls=":", label="mean")
ax1.set(xlabel="velocity [mm/s]", ylabel="radius [mm]", title="Hagen-Poiseuille profile")
ax1.legend()

# %%
# The fourth power of the radius
# ------------------------------
#
# Halving the bore cuts the flow sixteenfold, which is why narrowed arteries
# matter so much.

radii = np.geomspace(0.25e-3, 2e-3, 20)
ax2.loglog(1e3 * radii, poiseuille_flow(dp, radii, L, mu) * 1e9, marker="o")
ax2.set(xlabel="radius [mm]", ylabel="flow [mm³/s]", title="Q grows as R⁴")
fig.tight_layout()
slope = np.polyfit(np.log(radii), np.log(poiseuille_flow(dp, radii, L, mu)), 1)[0]
print(f"log-log slope of Q against R: {slope:.3f}")

# %%
# The same law in Darcy's form
# ----------------------------
#
# The profile integrates to the flow, and the Darcy head loss with
# :math:`f = 64/Re` reproduces the applied pressure drop.

Re = reynolds_number(rho, v_mean, 2 * R, mu)
head = darcy_head_loss(laminar_friction_factor(Re), L, 2 * R, v_mean)
print(f"Re = {Re:.1f}; Darcy loss {rho * G * head:.2f} Pa against {dp:.2f} Pa applied")

assert np.isclose(integrate.quad(lambda s: 2 * np.pi * s * poiseuille_velocity(s, dp, R, L, mu), 0, R)[0], Q)
assert np.isclose(u.max(), 2 * v_mean) and np.isclose(slope, 4.0)
assert np.isclose(rho * G * head, dp)

plt.show()
