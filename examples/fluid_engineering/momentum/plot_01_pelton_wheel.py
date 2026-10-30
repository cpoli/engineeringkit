r"""
Pelton's water wheel (1880)
===========================

Lester Pelton, a millwright in the California gold country, patented in
1880 a wheel whose buckets split the jet in two and turn each half almost
back on itself. By the momentum equation a vane that turns a jet through
:math:`\theta` is pushed with :math:`\rho A (V - u)^2 (1 - \cos\theta)`
(:func:`~engineeringkit.fluid_engineering.jet_force`); a wheel, which
always has a bucket in the jet, takes
:math:`\rho Q u (V - u)(1 - \cos\theta)`
(:func:`~engineeringkit.fluid_engineering.pelton_wheel_power`). Turning
the jet through nearly 180° and running the buckets at half the jet speed
leaves the water with almost no speed at all. This example compares
Pelton's buckets with the flat paddles of the hurdy-gurdy wheels they
replaced.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import G, jet_force, pelton_wheel_power

# %%
# A jet from 300 m of head
# ------------------------
#
# A 0.1 m nozzle with velocity coefficient 0.98.

H, d_jet, rho = 300.0, 0.1, 1000.0
V = 0.98 * np.sqrt(2 * G * H)
A = np.pi * d_jet**2 / 4
Q = A * V
jet_power = 0.5 * rho * Q * V**2
print(f"jet speed {V:.1f} m/s, flow {Q:.3f} m³/s, jet power {jet_power / 1e6:.2f} MW")

# %%
# Power against bucket speed
# --------------------------
#
# Flat paddles turn the jet through 90°; Pelton's buckets through 165°.

u = np.linspace(0, V, 200)
fig, ax = plt.subplots(figsize=(6, 4))
for label, theta in [("flat paddles, 90°", np.pi / 2), ("Pelton buckets, 165°", np.radians(165))]:
    eta = pelton_wheel_power(rho, Q, V, u, theta) / jet_power
    ax.plot(u / V, eta, label=label)
    print(f"{label}: best efficiency {eta.max():.3f} at u/V = {u[eta.argmax()] / V:.2f}")
ax.set(xlabel="bucket speed / jet speed", ylabel="hydraulic efficiency", title="Pelton wheel")
ax.legend()
fig.tight_layout()

# %%
# The force on a single vane
# --------------------------
#
# One bucket held still feels the full reversal force; one running away at
# half the jet speed feels a quarter of it, since it receives water at half
# the rate and turns it at half the relative speed.

theta = np.radians(165)
held, _ = jet_force(rho, A, V, theta)
running, _ = jet_force(rho, A, V, theta, vane_speed=V / 2)
print(f"bucket held still {held / 1e3:.1f} kN, running at V/2 {running / 1e3:.1f} kN")

assert np.isclose(pelton_wheel_power(rho, Q, V, V / 2, theta) / jet_power, (1 - np.cos(theta)) / 2)
assert np.isclose(running, held / 4)
assert pelton_wheel_power(rho, Q, V, V / 2, theta) > 1.9 * pelton_wheel_power(rho, Q, V, V / 2, np.pi / 2)

plt.show()
