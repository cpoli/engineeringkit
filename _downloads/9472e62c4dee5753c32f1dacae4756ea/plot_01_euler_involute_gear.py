r"""
Euler's involute gear tooth
===========================

Two gears turn at a constant speed ratio only if the common normal at
their point of contact always passes through one fixed point on the line of
centres, the *pitch point*. Leonhard Euler (1765) showed that the involute
of a circle, the curve traced by the end of a string unwound from it, meets
this condition. Its normal is always tangent to the base circle, so the
contact normal of two involutes is the fixed line tangent to both base
circles. The ratio is then the ratio of base radii, that is of tooth counts
(:func:`~engineeringkit.mechanical_engineering.gear_ratio`), and it does not
change if the centres are pulled apart. The example checks both properties
and the power balance of a gear train
(:func:`~engineeringkit.mechanical_engineering.gear_output_speed`,
:func:`~engineeringkit.mechanical_engineering.gear_output_torque`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import gear_output_speed, gear_output_torque, gear_ratio, power_torque

# %%
# The involute's normal is tangent to the base circle
# ---------------------------------------------------
#
# The involute of a base circle of radius :math:`r_b` is
# :math:`(r_b(\cos t + t\sin t),\ r_b(\sin t - t\cos t))`. At each point the
# normal line's distance from the gear centre should equal :math:`r_b`.

module, pressure_angle = 0.004, np.radians(20)
N1, N2 = 18, 45
rb1 = module * N1 / 2 * np.cos(pressure_angle)
rb2 = module * N2 / 2 * np.cos(pressure_angle)
t = np.linspace(0.05, 0.9, 200)
curve = rb1 * np.column_stack((np.cos(t) + t * np.sin(t), np.sin(t) - t * np.cos(t)))
tangent = np.gradient(curve, t, axis=0)
normal = np.column_stack((-tangent[:, 1], tangent[:, 0])) / np.hypot(*tangent.T)[:, None]
distance = np.abs(curve[:, 0] * normal[:, 1] - curve[:, 1] * normal[:, 0])
print(f"base radius {1e3 * rb1:.3f} mm; normal-line distance from centre {1e3 * distance.min():.3f}-{1e3 * distance.max():.3f} mm")
assert np.allclose(distance[1:-1], rb1, rtol=1e-4)

# %%
# The ratio survives a change of centre distance
# ----------------------------------------------
#
# The line of action is tangent to both base circles, so it cuts the line
# of centres :math:`C` at :math:`C\,r_{b1}/(r_{b1} + r_{b2})` from gear 1.
# The speed ratio, the inverse ratio of those pitch radii, stays
# :math:`N_2/N_1` while the working pressure angle
# :math:`\arccos((r_{b1} + r_{b2})/C)` grows.

nominal = module * (N1 + N2) / 2
for C in (nominal, 1.02 * nominal, 1.05 * nominal):
    pitch1 = C * rb1 / (rb1 + rb2)
    ratio = (C - pitch1) / pitch1
    angle = np.degrees(np.arccos((rb1 + rb2) / C))
    print(f"C = {1e3 * C:.2f} mm: speed ratio {ratio:.6f}, pressure angle {angle:.2f} deg")
    assert np.isclose(ratio, gear_ratio(N2, N1))

# %%
# Power through a two-stage train
# -------------------------------
#
# Speed falls and torque rises by the tooth ratio at each mesh; with 98%
# efficient meshes the output power is :math:`0.98^2` of the input.

speed, torque = 1450 * 2 * np.pi / 60, 12.0
p_in = power_torque(torque, speed)
for driver, driven in ((N1, N2), (20, 50)):
    speed, torque = gear_output_speed(speed, driver, driven), gear_output_torque(torque, driver, driven, efficiency=0.98)
print(f"output {speed * 60 / (2 * np.pi):.1f} rpm, {torque:.1f} N m, efficiency {power_torque(torque, speed) / p_in:.4f}")
assert np.isclose(power_torque(torque, speed) / p_in, 0.98**2)

fig, ax = plt.subplots(figsize=(5, 5))
circle = np.linspace(0, 2 * np.pi, 400)
ax.plot(rb1 * np.cos(circle), rb1 * np.sin(circle), "0.6", lw=1, label="base circle")
ax.plot(*curve.T, lw=2, label="involute flank")
for i in range(0, len(t), 25):
    foot = curve[i] - normal[i] * np.dot(curve[i], normal[i])
    ax.plot([curve[i, 0], foot[0]], [curve[i, 1], foot[1]], "C3", lw=0.8)
ax.plot([], [], "C3", lw=0.8, label="normals, tangent to base circle")
ax.set_aspect("equal")
ax.set(xlabel="x [m]", ylabel="y [m]", title=f"Involute of an {N1}-tooth gear")
ax.legend(loc="lower left", fontsize=8)
fig.tight_layout()

plt.show()
