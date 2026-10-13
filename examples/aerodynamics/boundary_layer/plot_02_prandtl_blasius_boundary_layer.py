r"""
Prandtl's boundary layer and the Blasius solution
=================================================

At the 1904 Heidelberg congress Ludwig Prandtl proposed that viscosity
matters only in a thin *boundary layer* next to a surface. His student Heinrich
Blasius (1908) reduced the flat-plate layer to the ordinary differential
equation :math:`f''' + \tfrac12 f f'' = 0`. Solving it by shooting
recovers the constants behind
:func:`~engineeringkit.aerodynamics.boundary_layer_laminar`
(:math:`\delta_{99} \approx 5x/\sqrt{Re_x}`) and
:func:`~engineeringkit.aerodynamics.skin_friction_laminar`
(:math:`C_f = 1.328/\sqrt{Re_L}`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from engineeringkit.aerodynamics import boundary_layer_laminar, reynolds_number, skin_friction_laminar

# %%
# Shoot for the wall shear f''(0)
# -------------------------------


def blasius(t, y):
    return [y[1], y[2], -0.5 * y[0] * y[2]]


def edge_velocity(curvature):
    return solve_ivp(blasius, (0, 12), [0, 0, curvature], rtol=1e-10, atol=1e-12).y[1, -1] - 1


fpp0 = brentq(edge_velocity, 0.1, 1.0)
eta = np.linspace(0, 8, 801)
sol = solve_ivp(blasius, (0, 8), [0, 0, fpp0], t_eval=eta, rtol=1e-10, atol=1e-12)
eta99 = np.interp(0.99, sol.y[1], eta)
print(f"f''(0) = {fpp0:.5f} -> local C_f sqrt(Re_x) = {2 * fpp0:.4f}, average = {4 * fpp0:.4f}")
print(f"delta_99 sqrt(Re_x) / x = {eta99:.3f}")
assert np.isclose(4 * fpp0, skin_friction_laminar(1.0), rtol=1e-3)
assert abs(eta99 - 4.91) < 0.01

# %%
# Growth of the layer on a 1 m plate
# ----------------------------------

x = np.linspace(1e-3, 1.0, 200)
Re_x = reynolds_number(1.225, 5.0, x, 1.789e-5)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
axes[0].plot(sol.y[1], eta)
axes[0].axhline(eta99, color="gray", ls="--", label=f"99% at eta = {eta99:.2f}")
axes[0].set(xlabel="u / U", ylabel="similarity variable eta", title="Blasius velocity profile")
axes[0].legend()
axes[1].plot(x, boundary_layer_laminar(x, Re_x) * 1e3, label="5 x / sqrt(Re_x)")
axes[1].plot(x, eta99 * x / np.sqrt(Re_x) * 1e3, "--", label="Blasius 4.91 x / sqrt(Re_x)")
axes[1].set(xlabel="distance from leading edge [m]", ylabel="thickness [mm]", title="Air at 5 m/s")
axes[1].legend()
fig.tight_layout()

plt.show()
