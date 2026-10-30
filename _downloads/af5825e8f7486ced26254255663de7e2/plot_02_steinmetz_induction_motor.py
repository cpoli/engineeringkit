r"""
Steinmetz's induction motor equivalent circuit (1897)
=====================================================

Tesla's induction motor had no electrical connection to its rotor, which
made it hard to calculate. Charles Proteus Steinmetz treated it as a
transformer whose secondary is shorted through a resistance
:math:`R_2/s` that depends on the slip :math:`s` between the rotating
field and the rotor. The power absorbed by :math:`R_2/s` is the air-gap
power, and the torque follows from it
(:func:`~engineeringkit.electrical_engineering.induction_motor_torque`).
The torque peaks at the breakdown slip, where :math:`R_2/s` matches the
impedance in front of it
(:func:`~engineeringkit.electrical_engineering.induction_motor_breakdown`),
and the peak itself does not depend on the rotor resistance. This example
draws the torque-speed curve of a four-pole, 50 Hz motor, checks the
breakdown point and shows why wound-rotor motors add resistance to start.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import induction_motor_breakdown, induction_motor_torque

# %%
# The torque-speed curve
# ----------------------
#
# A 400 V (231 V per phase), four-pole, 50 Hz motor: synchronous speed
# 1500 rpm.

v_phase, poles, f = 400 / np.sqrt(3), 4, 50.0
ws = 2 * np.pi * f / (poles / 2)
r1, x1, r2, x2, xm = 0.3, 0.6, 0.25, 0.6, 20.0
slip = np.linspace(1e-4, 1, 2000)
speed_rpm = (1 - slip) * ws * 60 / (2 * np.pi)
torque = induction_motor_torque(v_phase, slip, ws, r1, x1, r2, x2, xm)
s_max, t_max = induction_motor_breakdown(v_phase, ws, r1, x1, r2, x2, xm)
t_start = induction_motor_torque(v_phase, 1.0, ws, r1, x1, r2, x2, xm)
print(f"breakdown torque {t_max:.0f} N m at slip {s_max:.3f}; starting torque {t_start:.0f} N m")
assert np.isclose(t_max, torque.max(), rtol=1e-4)

# %%
# Rated operation
# ---------------
#
# At small slip the torque is nearly proportional to slip, so the speed
# droops only a few percent under load.

s_rated = 0.03
t_rated = induction_motor_torque(v_phase, s_rated, ws, r1, x1, r2, x2, xm)
p_out = t_rated * (1 - s_rated) * ws
print(f"at 3% slip: {t_rated:.0f} N m, {p_out / 1e3:.1f} kW mechanical (air-gap power minus rotor copper loss)")
assert np.isclose(p_out, t_rated * ws - s_rated * t_rated * ws)

# %%
# Starting with extra rotor resistance
# ------------------------------------
#
# Choosing :math:`R_2` so that the breakdown slip is 1 puts the peak torque
# at standstill.

r2_start = r2 / s_max
s_start, t_start_max = induction_motor_breakdown(v_phase, ws, r1, x1, r2_start, x2, xm)
print(f"rotor resistance {r2_start:.2f} ohm: breakdown at slip {s_start:.3f} with the same {t_start_max:.0f} N m")
assert np.isclose(s_start, 1) and np.isclose(t_start_max, t_max)

fig, ax = plt.subplots(figsize=(7, 3.8))
ax.plot(speed_rpm, torque, label=f"R2 = {r2} ohm")
ax.plot(speed_rpm, induction_motor_torque(v_phase, slip, ws, r1, x1, r2_start, x2, xm), label=f"R2 = {r2_start:.2f} ohm")
ax.plot((1 - s_max) * ws * 60 / (2 * np.pi), t_max, "ko")
ax.set(xlabel="speed [rpm]", ylabel="torque [N m]", title="Induction motor torque-speed curve")
ax.legend(fontsize="small")
fig.tight_layout()

plt.show()
