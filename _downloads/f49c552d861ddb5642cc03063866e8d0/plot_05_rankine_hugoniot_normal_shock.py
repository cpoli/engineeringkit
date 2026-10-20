r"""
Rankine-Hugoniot jump conditions across a normal shock
======================================================

W. J. Macquorn Rankine (1870) and Pierre-Henri Hugoniot (1887-1889) derived
the jump in state across a shock from conservation of mass, momentum and
energy alone. For a normal shock in a perfect gas the downstream Mach number
and pressure ratio are
:func:`~engineeringkit.aerodynamics.normal_shock_mach2` and
:func:`~engineeringkit.aerodynamics.normal_shock_pressure_ratio`. This example
reconstructs the full downstream state and checks the three conservation
laws, the entropy rise, and the loss of stagnation pressure.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import normal_shock_mach2, normal_shock_pressure_ratio, speed_of_sound, stagnation_pressure

gamma, R = 1.4, 287.05287
p1, T1 = 101325.0, 288.15
M1 = np.linspace(1.0, 5.0, 81)

# %%
# Reconstruct the downstream state
# --------------------------------

rho1 = p1 / (R * T1)
u1 = M1 * speed_of_sound(gamma, R, T1)
p2 = p1 * normal_shock_pressure_ratio(M1, gamma)
density_ratio = (gamma + 1) * M1**2 / ((gamma - 1) * M1**2 + 2)
rho2 = rho1 * density_ratio
T2 = p2 / (rho2 * R)
u2 = normal_shock_mach2(M1, gamma) * speed_of_sound(gamma, R, T2)
cp = gamma * R / (gamma - 1)

mass = rho2 * u2 / (rho1 * u1) - 1
momentum = (p2 + rho2 * u2**2) / (p1 + rho1 * u1**2) - 1
energy = (cp * T2 + u2**2 / 2) / (cp * T1 + u1**2 / 2) - 1
print(f"max conservation residuals: mass {np.max(abs(mass)):.1e}, momentum {np.max(abs(momentum)):.1e}, energy {np.max(abs(energy)):.1e}")
assert np.allclose([mass, momentum, energy], 0, atol=1e-12)

# %%
# Entropy rises and stagnation pressure falls
# -------------------------------------------

entropy = cp * np.log(T2 / T1) - R * np.log(p2 / p1)
p01 = stagnation_pressure(p1, M1, gamma)
p02 = stagnation_pressure(p2, normal_shock_mach2(M1, gamma), gamma)
assert np.all(entropy >= -1e-9)
print(f"at M1 = 2: p02/p01 = {np.interp(2, M1, p02 / p01):.4f}, M2 = {float(normal_shock_mach2(2.0)):.4f}")

fig, axes = plt.subplots(1, 3, figsize=(11, 3.3))
axes[0].plot(M1, normal_shock_mach2(M1, gamma))
axes[0].set(xlabel="$M_1$", ylabel="$M_2$", title="Downstream Mach")
axes[1].plot(M1, p2 / p1, label="$p_2/p_1$")
axes[1].plot(M1, density_ratio, label=r"$\rho_2/\rho_1$ (→ 6)")
axes[1].plot(M1, T2 / T1, label="$T_2/T_1$")
axes[1].set(xlabel="$M_1$", title="Static jumps")
axes[1].legend()
axes[2].plot(M1, p02 / p01)
axes[2].set(xlabel="$M_1$", ylabel="$p_{02}/p_{01}$", title="Stagnation-pressure loss")
fig.tight_layout()

plt.show()
