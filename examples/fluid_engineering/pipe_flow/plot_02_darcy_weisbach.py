r"""
Weisbach and Darcy: head loss in pipes
======================================

Julius Weisbach wrote the friction loss in a pipe in 1845 as
:math:`h_f = f\,(L/D)\,v^2/(2g)`, with a coefficient that he fitted to
experiments. Henry Darcy's tests on new and incrusted cast-iron mains in
Paris, published in 1857, showed that the coefficient depends on the
diameter and the state of the wall. The Darcy-Weisbach equation
(:func:`~engineeringkit.fluid_engineering.darcy_head_loss`), with loss
coefficients for the fittings
(:func:`~engineeringkit.fluid_engineering.minor_head_loss`), is still how
pipelines are sized. This example follows the energy and hydraulic grade
lines of water flowing from a reservoir through two pipes to a nozzle.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import (
    colebrook_friction_factor,
    darcy_head_loss,
    minor_head_loss,
    reynolds_number,
    velocity_head,
)
from engineeringkit.fluid_engineering.visualizers import plot_grade_lines

# %%
# The system
# ----------
#
# A reservoir surface stands 25 m above the outlet. 200 m of 150 mm pipe
# lead to a contraction and 100 m of 100 mm pipe, ending in a free jet.
# Water flows at 25 L/s; the entrance has :math:`K = 0.5` and the
# contraction :math:`K = 0.3` on the downstream velocity.

rho, mu, eps, Q, z_res = 998.0, 1.0e-3, 0.045e-3, 0.025, 25.0
pipes = [(200.0, 0.15), (100.0, 0.10)]
v = [Q / (np.pi * D**2 / 4) for _, D in pipes]
f = [colebrook_friction_factor(reynolds_number(rho, vi, D, mu), eps / D) for vi, (_, D) in zip(v, pipes)]

losses = {
    "entrance": minor_head_loss(0.5, v[0]),
    "pipe 1": darcy_head_loss(f[0], pipes[0][0], pipes[0][1], v[0]),
    "contraction": minor_head_loss(0.3, v[1]),
    "pipe 2": darcy_head_loss(f[1], pipes[1][0], pipes[1][1], v[1]),
}
for name, h in losses.items():
    print(f"{name:>11}: {h:5.2f} m")

# %%
# Grade lines
# -----------
#
# The energy grade line starts at the reservoir surface and drops at each
# loss; the hydraulic grade line sits one velocity head below it. The pipe
# runs along the ground at 0 m.

x = np.array([0, 0, 200, 200, 300])
egl = z_res - np.cumsum([losses["entrance"], 0, losses["pipe 1"], losses["contraction"], losses["pipe 2"]])
hv = velocity_head(np.array([v[0], v[0], v[0], v[1], v[1]]))
ax = plot_grade_lines(x, np.zeros(5), egl - hv, hv)
ax.set_title("Reservoir to free jet")
ax.figure.tight_layout()

# %%
# The energy left at the outlet is the jet's velocity head; the rest went
# to friction and fittings.

jet = velocity_head(v[1])
print(f"jet velocity head {jet:.2f} m + losses {sum(losses.values()):.2f} m = {jet + sum(losses.values()):.2f} m available")
print(f"friction share of the losses: {(losses['pipe 1'] + losses['pipe 2']) / sum(losses.values()):.0%}")

# Head loss in fully rough flow grows as v²: doubling the flow roughly quadruples it.
h2 = darcy_head_loss(colebrook_friction_factor(reynolds_number(rho, 2 * v[1], 0.1, mu), eps / 0.1), 100.0, 0.1, 2 * v[1])
print(f"doubling the flow multiplies the pipe-2 loss by {h2 / losses['pipe 2']:.2f}")

assert z_res - sum(losses.values()) > jet  # the pipe can carry this flow with head to spare
assert 3.6 < h2 / losses["pipe 2"] < 4.0
assert np.all(np.diff(egl) <= 0)

plt.show()
