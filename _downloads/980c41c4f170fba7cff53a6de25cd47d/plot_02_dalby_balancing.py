r"""
Dalby's balancing of rotating masses (1902)
===========================================

Faster engines and turbines made unbalance a practical danger: a mass
:math:`m` at radius :math:`r` pulls on the bearings with
:math:`mr\omega^2`, a force that grows with the square of speed. William
Dalby's *The Balancing of Engines* (1902) taught designers to balance any
set of rotating masses with two correction masses in chosen planes. Taking
moments about one correction plane, the *reference plane*, removes the
unknown correction there and gives the other directly; the force balance
then gives the first. A single correction can only make the shaft balance
statically, resting in any position, and leaves a rocking couple.
:func:`~engineeringkit.mechanical_engineering.balance_rotating_masses` finds
the corrections and :func:`~engineeringkit.mechanical_engineering.unbalance`
measures what remains. The example simulates the bearing reactions of a
rigid rotor before and after balancing.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import balance_rotating_masses, unbalance

# %%
# Four masses on a shaft
# ----------------------

masses = np.array([0.20, 0.30, 0.25, 0.15])
radii = np.array([0.08, 0.06, 0.07, 0.09])
angles = np.radians([0.0, 75.0, 160.0, 270.0])
z = np.array([0.15, 0.30, 0.45, 0.60])
planes = np.array([0.0, 0.75])
force, couple = unbalance(masses, radii, angles, z)
print(f"before: force unbalance {force * 1e3:.2f} g m, couple unbalance {couple * 1e3:.2f} g m^2")

static = balance_rotating_masses(masses, radii, angles)
dynamic = balance_rotating_masses(masses, radii, angles, z, planes)
print("static correction:", static.summary())
print("two-plane corrections:", dynamic.summary())

# %%
# Bearing reactions
# -----------------
#
# Bearings at the two correction planes carry the rotating forces of a
# rigid rotor. Static balance removes the net force but the two bearings
# still push in opposite directions; two-plane balance removes both.


def bearing_forces(m_list, r_list, a_list, z_list, omega):
    u = np.asarray(m_list) * r_list * np.exp(1j * np.asarray(a_list))
    right = (u * (np.asarray(z_list) - planes[0])).sum() / (planes[1] - planes[0]) * omega**2
    left = u.sum() * omega**2 - right
    return abs(left), abs(right)


omega = 2 * np.pi * 3000 / 60
before = bearing_forces(masses, radii, angles, z, omega)
once = bearing_forces([*masses, 1.0], [*radii, static.mass_radius[0]], [*angles, static.angle[0]], [*z, planes[0]], omega)
twice = bearing_forces([*masses, 1.0, 1.0], [*radii, *dynamic.mass_radius], [*angles, *dynamic.angle], [*z, *planes], omega)
for label, f in (("unbalanced", before), ("static balance", once), ("two-plane balance", twice)):
    print(f"{label:>17}: bearing forces {f[0]:8.2f} N and {f[1]:8.2f} N at 3000 rpm")
assert np.isclose(once[0], once[1]) and once[0] > 1
assert max(twice) < 1e-10
f_after, c_after = unbalance([*masses, 1.0, 1.0], [*radii, *dynamic.mass_radius], [*angles, *dynamic.angle], [*z, *planes])
assert f_after < 1e-15 and c_after < 1e-15

fig, ax = dynamic.plot()
for m_i, r_i, a_i in zip(masses, radii, angles):
    ax.plot([a_i, a_i], [0, m_i * r_i], "k-", lw=1)
ax.set_title("Masses (black) and two-plane corrections [kg m]")

plt.show()
