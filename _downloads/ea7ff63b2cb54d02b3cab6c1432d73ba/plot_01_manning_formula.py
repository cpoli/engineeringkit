r"""
Manning's open-channel formula
==============================

Robert Manning, chief engineer of the Irish Board of Works, compared the
open-channel formulas of his day against the best measurements in 1889 and
proposed that the mean velocity of uniform flow is
:math:`V = R^{2/3}S^{1/2}/n` in SI units
(:func:`~engineeringkit.fluid_engineering.manning_velocity`), with
:math:`R` the hydraulic radius
(:func:`~engineeringkit.fluid_engineering.hydraulic_radius`), :math:`S`
the bed slope and :math:`n` a roughness coefficient. It is Chézy's formula
:math:`V = C\sqrt{RS}` with :math:`C = R^{1/6}/n`. This example finds the
normal depth of a concrete channel and checks whether the flow is
subcritical.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from engineeringkit.fluid_engineering import froude_number, hydraulic_radius, manning_flow, manning_velocity

# %%
# Rating curve of a rectangular channel
# -------------------------------------
#
# A concrete channel (:math:`n = 0.013`) 3 m wide on a slope of 1 in 1000.

n, b, S = 0.013, 3.0, 1e-3
y = np.linspace(0.05, 2.5, 100)
Q = manning_flow(n, b * y, b + 2 * y, S)

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(Q, y)
ax.set(xlabel="flow [m³/s]", ylabel="normal depth [m]", title="Rating curve, n = 0.013, S = 0.001")
fig.tight_layout()

# %%
# Normal depth for a design flow
# ------------------------------
#
# Manning's formula is solved for the depth that carries 10 m³/s.

Q_design = 10.0
y_n = optimize.brentq(lambda d: manning_flow(n, b * d, b + 2 * d, S) - Q_design, 0.01, 10.0)
R = hydraulic_radius(b * y_n, b + 2 * y_n)
V = manning_velocity(n, R, S)
Fr = froude_number(V, y_n)
ax.plot(Q_design, y_n, "ko")
print(f"normal depth {y_n:.3f} m, R = {R:.3f} m, V = {V:.2f} m/s, Fr = {Fr:.2f}")
print("subcritical: disturbances travel upstream" if Fr < 1 else "supercritical")

# %%
# Chézy's coefficient
# -------------------
#
# Manning's :math:`C = R^{1/6}/n` changes slowly with depth, which is why
# Chézy's constant-C formula worked within a single river.

C = R ** (1 / 6) / n
print(f"Chézy C = {C:.1f} m^(1/2)/s")

assert np.isclose(manning_flow(n, b * y_n, b + 2 * y_n, S), Q_design)
assert np.isclose(V, C * np.sqrt(R * S))
assert Fr < 1 and np.all(np.diff(Q) > 0)

plt.show()
