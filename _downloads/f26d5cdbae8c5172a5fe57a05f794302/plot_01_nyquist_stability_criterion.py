r"""
Nyquist's stability criterion
=============================

Long-distance telephone repeaters with feedback sometimes "sang". Harry
Nyquist (1932) showed how to predict this from the measured open-loop
frequency response alone. Trace :math:`L(j\omega)` for all frequencies
(:func:`~engineeringkit.control_systems.nyquist`): the closed loop has
:math:`Z = N + P` unstable poles, where :math:`N` counts clockwise
encirclements of :math:`-1`
(:func:`~engineeringkit.control_systems.nyquist_encirclements`) and
:math:`P` the unstable open-loop poles.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems import feedback, nyquist_encirclements, poles, transfer_function
from engineeringkit.control_systems.visualizers import plot_nyquist

# %%
# A stable plant at two gains
# ---------------------------
#
# :math:`K/(s+1)^3` crosses the negative real axis at :math:`-K/8`.

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for K in [4.0, 12.0]:
    loop = transfer_function([K], [1, 3, 3, 1])
    N = nyquist_encirclements(loop)
    Z = int(np.sum(poles(feedback(loop)).real > 0))
    print(f"K = {K:4.1f}: N = {N}, P = 0, closed-loop unstable poles Z = {Z}")
    assert Z == N + 0
    plot_nyquist(loop, ax=axes[0], label=f"K = {K:g}, N = {N}")
axes[0].set(xlim=(-3, 13), ylim=(-8, 8), title=r"$K/(s+1)^3$")

# %%
# An unstable plant stabilised by feedback
# ----------------------------------------
#
# :math:`2/(s-1)` has one unstable pole (:math:`P = 1`); one
# counter-clockwise encirclement (:math:`N = -1`) gives :math:`Z = 0`.

loop = transfer_function([2], [1, -1])
N = nyquist_encirclements(loop)
Z = int(np.sum(poles(feedback(loop)).real > 0))
print(f"2/(s-1): N = {N}, P = 1, Z = {Z}")
assert N == -1 and Z == N + 1
plot_nyquist(loop, ax=axes[1], label="$2/(s-1)$")
axes[1].set(title="Unstable open loop, stable closed loop")
axes[1].set_aspect("equal")
fig.tight_layout()

plt.show()
