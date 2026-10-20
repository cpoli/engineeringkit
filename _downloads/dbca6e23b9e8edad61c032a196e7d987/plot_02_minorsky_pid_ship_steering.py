r"""
Minorsky's PID controller for ship steering
===========================================

Nicolas Minorsky (1922) watched helmsmen on the battleship *New Mexico* and
noticed they reacted to the heading error, to how fast it was changing, and
to any error that persisted. His three-term controller,
:math:`\delta = K_p e + K_i \int e\,dt + K_d \dot e`, is the
:class:`~engineeringkit.control_systems.PID`. Using Nomoto's ship model
:math:`\psi/\delta = K/(s(Ts + 1))` and a steady wind moment, proportional
plus derivative action leaves a constant heading offset :math:`d/K_p`;
adding integral action removes it.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems import PID, feedback, poles, step_response, transfer_function

K, T = 0.5, 10.0  # rudder gain [1/s], yaw time constant [s]
ship = transfer_function([K], [T, 1, 0])
wind = 0.05  # equivalent steady rudder offset [rad]

# %%
# Heading error caused by a steady wind
# -------------------------------------

t = np.linspace(0, 400, 4001)
controllers = {"PD": PID(2.0, 0.0, 20.0), "PID": PID(2.0, 0.1, 20.0)}
fig, ax = plt.subplots()
for name, pid in controllers.items():
    disturbance_to_heading = feedback(ship, pid.transfer_function())
    assert np.all(poles(disturbance_to_heading).real < 0)
    _, psi = step_response(disturbance_to_heading, t)
    ax.plot(t, np.degrees(wind * psi), label=name)
    print(f"{name}: final heading error {np.degrees(wind * psi[-1]):.3f} deg")
    expected = wind / pid.kp if pid.ki == 0 else 0.0
    assert np.isclose(wind * psi[-1], expected, atol=2e-4)
ax.axhline(np.degrees(wind / 2.0), color="gray", ls=":", label="$d/K_p$")
ax.set(xlabel="time [s]", ylabel="heading error [deg]", title="Nomoto ship model with a steady wind moment")
ax.legend()

plt.show()
