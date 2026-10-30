r"""
The Prandtl-Meyer expansion fan (1907-1908)
===========================================

Ludwig Prandtl (1907) and his student Theodor Meyer (1908) showed that a
supersonic stream turning *away* from itself around a corner expands
isentropically through a centred fan of Mach waves. The turning angle that
takes a sonic stream to Mach :math:`M` is the Prandtl-Meyer function
:math:`\nu(M)` (:func:`~engineeringkit.aerodynamics.prandtl_meyer_angle`),
so a turn by :math:`\theta` sets :math:`\nu(M_2) = \nu(M_1) + \theta`
(:func:`~engineeringkit.aerodynamics.prandtl_meyer_expansion`). Combined
with Meyer's oblique shocks this gives *shock-expansion theory*, which this
example applies to a flat plate at incidence and checks against Ackeret's
linear supersonic lift :math:`c_l = 4\alpha/\sqrt{M^2-1}`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import (
    oblique_shock,
    prandtl_meyer_angle,
    prandtl_meyer_expansion,
    stagnation_pressure,
)

gamma = 1.4

# %%
# Expansion round a 20° corner at Mach 1.5
# ----------------------------------------

M2 = float(prandtl_meyer_expansion(1.5, np.radians(20.0)))
p_ratio = float(stagnation_pressure(1.0, 1.5) / stagnation_pressure(1.0, M2))
print(f"nu(1.5) = {np.degrees(prandtl_meyer_angle(1.5)):.2f} deg -> M2 = {M2:.3f}, p2/p1 = {p_ratio:.4f}")
assert np.isclose(np.degrees(prandtl_meyer_angle(M2) - prandtl_meyer_angle(1.5)), 20.0)

# %%
# Shock-expansion theory for a flat plate
# ---------------------------------------
#
# The upper surface turns away from the stream through a fan and the lower
# surface turns into it through an oblique shock. The pressure difference
# acts normal to the plate.


def flat_plate_cl(M, alpha):
    p_inf = 1.0
    p0 = stagnation_pressure(p_inf, M, gamma)
    upper = stagnation_pressure(1.0, prandtl_meyer_expansion(M, alpha, gamma), gamma)
    p_upper = p0 / upper
    p_lower = p_inf * oblique_shock(M, alpha, gamma).pressure_ratio
    q = gamma / 2 * p_inf * M**2
    return (p_lower - p_upper) / q * np.cos(alpha)


alpha = np.radians(np.linspace(0.5, 15, 30))
fig, ax = plt.subplots(figsize=(7, 4))
for M in (2.0, 3.0, 4.0):
    exact = np.array([flat_plate_cl(M, a) for a in alpha])
    (line,) = ax.plot(np.degrees(alpha), exact, label=f"shock-expansion, M = {M:g}")
    ax.plot(np.degrees(alpha), 4 * alpha / np.sqrt(M**2 - 1), "--", color=line.get_color())
ax.set(xlabel="angle of attack [deg]", ylabel="$c_l$", title="Flat plate: shock-expansion (solid) vs linear theory (dashed)")
ax.legend()
fig.tight_layout()

small = np.radians(1.0)
print(f"M = 2, alpha = 1 deg: c_l = {flat_plate_cl(2.0, small):.5f}, linear theory {4 * small / np.sqrt(3):.5f}")
assert np.isclose(flat_plate_cl(2.0, small), 4 * small / np.sqrt(3), rtol=0.01)

plt.show()
