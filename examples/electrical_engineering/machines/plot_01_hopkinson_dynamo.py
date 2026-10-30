r"""
The Hopkinsons' dynamo and the DC machine (1886)
================================================

John Hopkinson introduced the *characteristic curve* of a dynamo, its EMF
against field current at constant speed, in 1879, and with his brother
Edward showed in 1886 how to predict it from the magnetic circuit of the
iron. Underneath both lies the armature equation of every DC machine: an
EMF :math:`E = k\Phi\omega` proportional to flux and speed, behind the
armature resistance :math:`R_a`, with torque :math:`T = k\Phi I_a`. As a
motor on a fixed supply the torque-speed curve is then a straight line from
the stall torque to the no-load speed
(:func:`~engineeringkit.electrical_engineering.dc_motor_torque`,
:func:`~engineeringkit.electrical_engineering.dc_motor_speed`); driven
faster than that the same machine generates. This example draws the
torque-speed lines for three armature voltages and checks the power
balance at an operating point.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import dc_motor_speed, dc_motor_torque

# %%
# Torque-speed lines
# ------------------
#
# A 240 V, 0.4 ohm armature with :math:`k\Phi = 2.2` V s/rad.

Ra, kphi = 0.4, 2.2
omega = np.linspace(0, 130, 200)
fig, ax = plt.subplots(figsize=(7, 3.8))
for V in (120.0, 180.0, 240.0):
    T = dc_motor_torque(V, Ra, kphi, omega)
    ax.plot(omega, T, label=f"V = {V:.0f} V")
    print(f"V = {V:.0f} V: no-load speed {V / kphi:.1f} rad/s, stall torque {kphi * V / Ra:.0f} N m")
assert dc_motor_torque(240.0, Ra, kphi, 240.0 / kphi) == 0

# %%
# Power balance at a load
# -----------------------
#
# Electrical input = copper loss + mechanical output, since
# :math:`VI_a = I_a^2R_a + EI_a` and :math:`EI_a = T\omega`.

V, T_load = 240.0, 150.0
w = dc_motor_speed(V, Ra, kphi, T_load)
Ia = T_load / kphi
p_in, p_copper, p_out = V * Ia, Ia**2 * Ra, T_load * w
print(f"at {T_load:.0f} N m: speed {w:.2f} rad/s ({w * 60 / (2 * np.pi):.0f} rpm), efficiency {p_out / p_in:.1%}")
assert np.isclose(p_in, p_copper + p_out)

# %%
# Generating
# ----------
#
# Above the no-load speed the EMF exceeds the supply, the current and torque
# reverse, and the machine returns power to the supply.

assert dc_motor_torque(V, Ra, kphi, 1.1 * V / kphi) < 0

ax.axhline(0, color="gray", lw=0.6)
ax.plot(w, T_load, "ko")
ax.set(xlabel="speed [rad/s]", ylabel="torque [N m]", title="Separately excited DC machine")
ax.legend(fontsize="small")
fig.tight_layout()

plt.show()
