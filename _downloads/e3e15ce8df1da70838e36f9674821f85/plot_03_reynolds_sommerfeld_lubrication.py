r"""
Reynolds and Sommerfeld: hydrodynamic lubrication
=================================================

Beauchamp Tower (1883) drilled an oil hole in the cap of a test bearing and
found the oil pushing the plug out: the film carried the load on a pressure
far above atmospheric. Osborne Reynolds (1886) explained it. A journal
running off-centre drags oil into a converging wedge, and in a thin film
the pressure obeys

.. math:: \frac{d}{d\theta}\left(h^3\frac{dp}{d\theta}\right) = 6\mu\omega r^2\frac{dh}{d\theta},
          \qquad h = c(1 + \varepsilon\cos\theta).

Arnold Sommerfeld (1904) integrated it in closed form for a long bearing.
The load, viscosity, speed and geometry then enter only through the
Sommerfeld number :math:`S = (r/c)^2\mu N/P`
(:func:`~engineeringkit.mechanical_engineering.sommerfeld_number`), which
fixes the eccentricity
(:func:`~engineeringkit.mechanical_engineering.long_bearing_eccentricity`).
The example integrates Reynolds's equation numerically and recovers
Sommerfeld's load and friction formulas.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

from engineeringkit.mechanical_engineering import long_bearing_eccentricity, petroff_friction_torque, sommerfeld_number

# %%
# Integrate Reynolds's equation round the journal
# -----------------------------------------------
#
# Integrating once gives :math:`dp/d\theta = 6\mu\omega r^2 (h - h^*)/h^3`.
# The constant :math:`h^*` makes the pressure return to its starting value
# after one turn. The film force and the shear on the journal follow by
# integrating round the circumference, per metre of bearing length.

viscosity, omega, r, L, c = 0.04, 2 * np.pi * 40, 0.025, 1.0, 30e-6
theta = np.linspace(0, 2 * np.pi, 20001)


def film(eps):
    h = c * (1 + eps * np.cos(theta))
    h_star = trapezoid(h**-2, theta) / trapezoid(h**-3, theta)
    dp = 6 * viscosity * omega * r**2 * (h - h_star) / h**3
    p = cumulative_trapezoid(dp, theta, initial=0)
    p -= p.mean()
    load = np.hypot(trapezoid(p * np.cos(theta), theta), trapezoid(p * np.sin(theta), theta)) * r * L
    torque = trapezoid(viscosity * omega * r / h + h / (2 * r) * dp, theta) * r**2 * L
    return p, load, torque


for eps in (0.2, 0.5, 0.8):
    p, load, torque = film(eps)
    W_sommerfeld = 12 * np.pi * viscosity * omega * r**3 * L * eps / (c**2 * (2 + eps**2) * np.sqrt(1 - eps**2))
    T_sommerfeld = 4 * np.pi * viscosity * omega * r**3 * L * (1 + 2 * eps**2) / (c * (2 + eps**2) * np.sqrt(1 - eps**2))
    S = sommerfeld_number(viscosity, omega, load, r, L, c)
    print(f"eps = {eps}: load {load / 1e3:7.1f} kN per metre (Sommerfeld {W_sommerfeld / 1e3:7.1f})")
    print(f"          S = {S:.4f} -> eps = {long_bearing_eccentricity(S):.4f}")
    assert np.isclose(load, W_sommerfeld, rtol=1e-6) and np.isclose(torque, T_sommerfeld, rtol=1e-6)
    assert np.isclose(long_bearing_eccentricity(S), eps)

# %%
# Lightly loaded bearings tend to Petroff's
# -----------------------------------------
#
# As the Sommerfeld number grows the journal moves to the centre and the
# friction torque falls to Petroff's value.

S = np.geomspace(0.01, 10, 60)
eps = long_bearing_eccentricity(S)
torque_ratio = (1 + 2 * eps**2) * 2 / ((2 + eps**2) * np.sqrt(1 - eps**2))
_, _, torque = film(0.01)
print(f"eps = 0.01: torque / Petroff = {torque / petroff_friction_torque(viscosity, r, L, c, omega):.6f}")
assert np.isclose(torque, petroff_friction_torque(viscosity, r, L, c, omega), rtol=1e-3)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for e in (0.2, 0.5, 0.8):
    axes[0].plot(np.degrees(theta), film(e)[0] / 1e6, label=f"eps = {e}")
axes[0].set(xlabel="angle from widest gap [deg]", ylabel="film pressure [MPa]", title="Full Sommerfeld pressure")
axes[0].legend()
axes[1].semilogx(S, eps, label="eccentricity ratio")
axes[1].semilogx(S, torque_ratio, label="torque / Petroff torque")
axes[1].set(xlabel="Sommerfeld number S", title="Long journal bearing")
axes[1].legend()
fig.tight_layout()

plt.show()
