r"""
James Thomson's triangular notch (1858)
=======================================

James Thomson, engineer and elder brother of Lord Kelvin, reported to the
British Association in 1858 that a V-shaped notch measures small flows
better than a rectangular weir. The flow section keeps the same shape as
the head changes, so the discharge coefficient hardly varies, and the flow
grows as :math:`H^{5/2}`
(:func:`~engineeringkit.fluid_engineering.v_notch_weir_flow`) instead of
:math:`H^{3/2}` (:func:`~engineeringkit.fluid_engineering.rectangular_weir_flow`),
so a small flow still gives a readable head. This example compares the
rating curves of the two weirs and finds how precisely each must be read.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import rectangular_weir_flow, rehbock_discharge_coefficient, v_notch_weir_flow

# %%
# Rating curves
# -------------
#
# A 90° notch and a 0.5 m rectangular weir standing 0.4 m above the bed of
# a laboratory flume, with Rehbock's coefficient for the rectangle.

H = np.linspace(0.005, 0.3, 200)
P, b = 0.4, 0.5
Q_vee = v_notch_weir_flow(H)
Q_rect = rectangular_weir_flow(H, b, rehbock_discharge_coefficient(H, P))

fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(Q_vee * 1000, H, label="90° V-notch, $Q \\propto H^{5/2}$")
ax.loglog(Q_rect * 1000, H, label="0.5 m rectangular weir, $Q \\propto H^{3/2}$")
ax.set(xlabel="flow [L/s]", ylabel="head over the crest [m]", title="Sharp-crested weir rating curves")
ax.legend()
fig.tight_layout()

# %%
# Head for a small flow
# ---------------------
#
# For 2 L/s the notch stands several times higher than the rectangle, so a
# 1 mm reading error matters far less.

Q = 2e-3
h_vee = np.interp(Q, Q_vee, H)
h_rect = np.interp(Q, Q_rect, H)
for name, h, n in [("V-notch", h_vee, 2.5), ("rectangular", h_rect, 1.5)]:
    print(f"{name}: head {h * 1000:.1f} mm; a 1 mm error is a {100 * n * 1e-3 / h:.1f}% flow error")

# %%
# The slopes of the log-log curves are the exponents of the head.

slope_vee = np.polyfit(np.log(H), np.log(Q_vee), 1)[0]
print(f"fitted exponent of the notch: {slope_vee:.3f}")

assert np.isclose(slope_vee, 2.5)
assert h_vee > 2.5 * h_rect

plt.show()
