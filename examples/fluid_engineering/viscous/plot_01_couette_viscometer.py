r"""
Couette's rotating-cylinder viscometer (1890)
=============================================

Maurice Couette, in his 1890 thesis at Paris, measured the viscosity of
water in the narrow gap between a fixed inner cylinder and a rotating outer
one, from the torque on the inner cylinder. Rotating the outer cylinder
keeps the flow laminar to high speeds. For a narrow gap the flow is the
plane shear flow between a fixed and a sliding plate, with a straight
velocity profile and uniform stress :math:`\mu U/h`
(:func:`~engineeringkit.fluid_engineering.couette_velocity`,
:func:`~engineeringkit.fluid_engineering.couette_shear_stress`); the exact
torque for any gap is
:func:`~engineeringkit.fluid_engineering.couette_viscometer_torque`. This
example reads a viscosity back from a torque and adds a pressure gradient
to the plane flow.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import couette_flow_rate, couette_shear_stress, couette_velocity, couette_viscometer_torque

# %%
# Reading a viscosity
# -------------------
#
# A bob of radius 20 mm and wetted length 60 mm in a cup of radius 21 mm,
# the cup turning at 60 rpm, reads a torque of 4.2 mN m on a glycerine
# solution. The torque is linear in :math:`\mu`, so one reading gives it.

R1, R2, L, omega = 0.020, 0.021, 0.060, 2 * np.pi
torque = 4.2e-3
mu = torque / couette_viscometer_torque(1.0, omega, R1, R2, L)
narrow = couette_shear_stress(mu, omega * R2, R2 - R1) * 2 * np.pi * R1 * L * R1
print(f"viscosity {mu:.4f} Pa s; the narrow-gap estimate of the torque is {narrow * 1e3:.2f} mN m")

# %%
# Couette-Poiseuille flow
# -----------------------
#
# Between flat plates 1 mm apart, the top one sliding at 1 m/s, a pressure
# gradient bends the line into a parabola. An adverse gradient beyond
# :math:`2\mu U/h^2` makes the fluid run backward near the fixed wall.

h, U = 1e-3, 1.0
y = np.linspace(0, h, 100)
critical = 2 * mu * U / h**2
fig, ax = plt.subplots(figsize=(5, 4))
for k in [-1.0, 0.0, 1.0, 3.0]:
    ax.plot(couette_velocity(y, h, U, mu, k * critical), y * 1000, label=f"dp/dx = {k:g} × 2μU/h²")
ax.axvline(0, color="k", lw=0.5)
ax.set(xlabel="velocity [m/s]", ylabel="height above fixed plate [mm]", title="Couette flow")
ax.legend(fontsize="small")
fig.tight_layout()
print(f"flow per metre width with no gradient: {couette_flow_rate(h, U, mu) * 1000:.2f} L/s")

assert np.isclose(couette_viscometer_torque(mu, omega, R1, R2, L), torque)
assert abs(narrow / torque - 1) < 0.05
assert np.isclose(couette_flow_rate(h, U, mu, critical * 3), U * h / 2 - h**3 * 3 * critical / (12 * mu))
assert couette_velocity(h / 10, h, U, mu, 3 * critical) < 0

plt.show()
