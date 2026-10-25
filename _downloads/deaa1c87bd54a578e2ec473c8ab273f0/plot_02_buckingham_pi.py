r"""
Buckingham's Pi theorem (1914)
==============================

Edgar Buckingham, at the US Bureau of Standards, published in 1914 the
theorem that turned dimensional analysis into a method: a physical
relation among :math:`n` variables whose dimensions involve :math:`r`
independent base units can be rewritten as a relation among :math:`n - r`
dimensionless products. The products are the null space of the dimension
matrix, which
:func:`~engineeringkit.fluid_engineering.buckingham_pi` finds by exact
row reduction. This example derives the groups of pipe friction and of
droplet break-up, then shows that pipes of different size and fluid with
the same groups have the same Euler number
(:func:`~engineeringkit.fluid_engineering.euler_number`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import (
    buckingham_pi,
    colebrook_friction_factor,
    darcy_head_loss,
    euler_number,
    reynolds_number,
    weber_number,
)

# %%
# Pressure drop in a pipe
# -----------------------
#
# Seven variables (pressure drop, density, speed, diameter, viscosity,
# roughness, length) in three base units give four groups. With density,
# speed and diameter repeating, they are the Euler number, the inverse
# Reynolds number, the relative roughness and the length ratio.

pipe = buckingham_pi(
    {
        "dp": {"M": 1, "L": -1, "T": -2},
        "rho": {"M": 1, "L": -3},
        "V": {"L": 1, "T": -1},
        "D": {"L": 1},
        "mu": {"M": 1, "L": -1, "T": -1},
        "eps": {"L": 1},
        "L": {"L": 1},
    },
    repeating=["rho", "V", "D"],
)
print(pipe.summary())
fig, ax = pipe.plot()
fig.set_size_inches(6, 3.2)
fig.tight_layout()

# %%
# Similar pipes, equal Euler numbers
# ----------------------------------
#
# A water main, an oil line and a small air duct, each at
# :math:`Re = 10^5`, :math:`\varepsilon/D = 10^{-4}` and :math:`L/D = 500`.
# The pressure drops differ by orders of magnitude; the Euler numbers agree.

cases = {"water": (998.0, 1.0e-3, 0.3), "oil": (870.0, 0.03, 0.5), "air": (1.2, 1.8e-5, 0.05)}
eu = []
for name, (rho, mu, D) in cases.items():
    V = 1e5 * mu / (rho * D)
    Re = reynolds_number(rho, V, D, mu)
    dp = rho * 9.80665 * darcy_head_loss(colebrook_friction_factor(Re, 1e-4), 500 * D, D, V)
    eu.append(euler_number(dp, rho, V))
    print(f"{name:>5}: V = {V:6.3f} m/s, dp = {dp:10.1f} Pa, Eu = {eu[-1]:.4f}")

# %%
# Droplet break-up
# ----------------
#
# A drop of diameter :math:`d` in a gas stream: drop and gas densities,
# speed, surface tension and the drop viscosity. The groups include the
# Weber number, which must exceed about 12 for the drop to shatter.

drop = buckingham_pi(
    {
        "rho_g": {"M": 1, "L": -3},
        "V": {"L": 1, "T": -1},
        "d": {"L": 1},
        "sigma": {"M": 1, "T": -2},
        "rho_l": {"M": 1, "L": -3},
        "mu_l": {"M": 1, "L": -1, "T": -1},
    }
)
print(drop.groups())
We = weber_number(1.2, 40.0, 1e-3, 0.072)
print(f"a 1 mm water drop in air at 40 m/s: We = {We:.1f}, {'breaks up' if We > 12 else 'survives'}")

assert not np.any(pipe.dimension_matrix @ pipe.exponents.T)
assert pipe.summary()["groups"] == 4 and len(drop.groups()) == 3
assert np.allclose(eu, eu[0])
assert "rho_g^-1 V^-2 d^-1 sigma^1" in drop.groups()

plt.show()
