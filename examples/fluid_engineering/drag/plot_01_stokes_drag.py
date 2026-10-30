r"""
Stokes's drag on a sphere (1851)
================================

George Gabriel Stokes, correcting pendulum experiments for the air dragged
along with the bob, solved the slow viscous flow past a sphere in a paper
read to the Cambridge Philosophical Society in December 1850 and printed in
1851. The drag is :math:`3\pi\mu dV`
(:func:`~engineeringkit.fluid_engineering.stokes_drag`), so a small sphere
settles at :math:`(\rho_p - \rho)gd^2/(18\mu)`
(:func:`~engineeringkit.fluid_engineering.stokes_settling_velocity`). The
law fails once inertia matters, near :math:`Re = 1`; the measured drag
curve (:func:`~engineeringkit.fluid_engineering.sphere_drag_coefficient`)
takes over. This example draws the curve and finds how fast sand grains of
different sizes fall through water
(:func:`~engineeringkit.fluid_engineering.terminal_velocity`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import (
    reynolds_number,
    sphere_drag_coefficient,
    stokes_drag,
    stokes_settling_velocity,
    terminal_velocity,
)

# %%
# The sphere drag curve
# ---------------------
#
# Stokes's :math:`C_D = 24/Re` is the straight line at the left; the curve
# levels off near 0.4 in Newton's range before the drag crisis.

Re = np.geomspace(0.01, 2e5, 200)
fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(Re, sphere_drag_coefficient(Re), label="drag curve (White's fit)")
ax.loglog(Re[Re < 10], 24 / Re[Re < 10], "--", label="Stokes, 24/Re")
ax.set(xlabel="Reynolds number", ylabel="drag coefficient $C_D$", title="Drag on a smooth sphere")
ax.legend()
ax.grid(True, which="both", lw=0.3)
fig.tight_layout()

# %%
# Settling sand
# -------------
#
# Quartz sand (2650 kg/m³) falling through water at 20 °C. Stokes's
# velocity is accurate for silt but overestimates coarse sand badly.

rho, mu, rho_p = 998.0, 1.0e-3, 2650.0
for d in [20e-6, 100e-6, 500e-6, 2e-3]:
    V = terminal_velocity(d, rho_p, rho, mu)
    V_stokes = stokes_settling_velocity(d, rho_p, rho, mu)
    print(f"d = {d * 1e6:6.0f} µm: V = {V * 1000:7.3f} mm/s, Stokes {V_stokes * 1000:9.3f} mm/s, Re = {reynolds_number(rho, V, d, mu):.3g}")

# %%
# The drag on a 20 µm silt grain at its settling speed equals its weight
# less its buoyancy.

d = 20e-6
V = stokes_settling_velocity(d, rho_p, rho, mu)
weight = (rho_p - rho) * 9.80665 * np.pi * d**3 / 6

assert np.isclose(stokes_drag(mu, V, d), weight)
assert np.isclose(terminal_velocity(d, rho_p, rho, mu), V, rtol=1e-2)
assert terminal_velocity(2e-3, rho_p, rho, mu) < 0.2 * stokes_settling_velocity(2e-3, rho_p, rho, mu)

plt.show()
