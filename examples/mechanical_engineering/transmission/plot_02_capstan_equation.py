r"""
The Euler-Eytelwein belt-friction law (1769)
============================================

A sailor holds a heavy load with a few turns of rope round a bollard.
Leonhard Euler (1762, published 1769) explained why: on each small arc
:math:`d\theta` of contact the rope presses on the drum with force
:math:`T\,d\theta`, friction can resist up to :math:`\mu T\,d\theta`, so
tension can grow by :math:`dT = \mu T\,d\theta` and
:math:`T_1/T_2 = e^{\mu\theta}`. Johann Albert Eytelwein (1808) brought the
result to engineers, for whom it sets the power a belt drive can transmit
(:func:`~engineeringkit.mechanical_engineering.belt_tension_ratio`,
:func:`~engineeringkit.mechanical_engineering.belt_speed`). The example
builds the exponential from small arcs, then sizes a flat and a V-belt.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import belt_speed, belt_tension_ratio

# %%
# Small arcs add up to an exponential
# -----------------------------------
#
# Each arc multiplies the tension by :math:`1 + \mu\,\Delta\theta`.

mu, wrap = 0.3, np.pi
for n in (4, 40, 400, 40000):
    print(f"{n:6d} arcs: T1/T2 = {(1 + mu * wrap / n) ** n:.5f}")
print(f"  exact     : T1/T2 = {belt_tension_ratio(mu, wrap):.5f}")
assert np.isclose((1 + mu * wrap / 40000) ** 40000, belt_tension_ratio(mu, wrap), rtol=1e-4)

# %%
# Turns round a bollard
# ---------------------

for turns in (0.5, 1, 2, 3):
    print(f"{turns:3} turns: a 100 N pull holds {100 * belt_tension_ratio(mu, 2 * np.pi * turns) / 1e3:7.2f} kN")

# %%
# Power a belt drive can carry
# ----------------------------
#
# A 200 mm pulley at 1450 rpm with 170 degrees of wrap and a tight side
# limited to 1.2 kN transmits :math:`P = (T_1 - T_2) v` before slipping. A
# V-belt in a 38-degree groove wedges in and grips :math:`1/\sin 19^\circ`
# times harder.

v = belt_speed(0.2, 1450)
T1, wrap = 1200.0, np.radians(170)
for name, ratio in (("flat belt", belt_tension_ratio(mu, wrap)), ("V-belt", belt_tension_ratio(mu, wrap, groove_angle=np.radians(38)))):
    print(f"{name:9s}: T1/T2 = {ratio:6.2f}, power {(T1 - T1 / ratio) * v / 1e3:.2f} kW at {v:.2f} m/s")
assert np.isclose(np.log(belt_tension_ratio(mu, wrap, groove_angle=np.radians(38))), mu * wrap / np.sin(np.radians(19)))

angle = np.linspace(0, 4 * np.pi, 400)
fig, ax = plt.subplots(figsize=(6, 3.5))
for coefficient in (0.1, 0.2, 0.3):
    ax.semilogy(angle / (2 * np.pi), belt_tension_ratio(coefficient, angle), label=f"mu = {coefficient}")
ax.set(xlabel="turns of wrap", ylabel="T1 / T2", title="Capstan amplification")
ax.legend()
fig.tight_layout()

plt.show()
