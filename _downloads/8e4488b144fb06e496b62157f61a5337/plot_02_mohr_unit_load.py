r"""
Mohr's unit-load method (1874)
==============================

To find how far one joint of a truss moves, Otto Mohr (1874) applied the
principle of virtual work with an imaginary unit load at that joint. The
internal forces :math:`n` it produces, working through the real member
extensions :math:`NL/(AE)`, must equal the work of the unit load through
the wanted displacement:

.. math::

    1\cdot\delta = \sum \frac{nNL}{AE}

(:func:`~engineeringkit.structural_engineering.unit_load_truss_deflection`).
Both sets of forces come from statics alone, so a determinate truss needs
no stiffness matrix. For beams the sum becomes :math:`\int mM/(EI)\,dx`
(:func:`~engineeringkit.structural_engineering.unit_load_beam_deflection`).
The example finds the midspan deflection of a Warren truss member by
member and checks it against the direct stiffness method.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import PlanarTruss, method_of_joints, unit_load_beam_deflection, unit_load_truss_deflection

# %%
# A Warren truss carrying two panel loads
# ---------------------------------------

nodes = np.array([[0, 0], [4, 0], [8, 0], [12, 0], [2, 3], [6, 3], [10, 3]], float)
members = ((0, 1), (1, 2), (2, 3), (0, 4), (4, 1), (1, 5), (5, 2), (2, 6), (6, 3), (4, 5), (5, 6))
fixed = [0, 1, 7]
E, A = 200e9, 2e-3
loads = np.zeros(14)
loads[[3, 5]] = -60e3
unit = np.zeros(14)
unit[5] = -1.0  # downward unit load at the midspan joint 2

N = method_of_joints(nodes, members, loads, fixed).member_forces
n = method_of_joints(nodes, members, unit, fixed).member_forces
lengths = np.array([np.hypot(*(nodes[j] - nodes[i])) for i, j in members])
contributions = n * N * lengths / (A * E)
delta = unit_load_truss_deflection(N, n, lengths, A, E)
print(f"midspan deflection by virtual work: {delta * 1e3:.3f} mm")

stiffness = PlanarTruss(nodes, members, E, A).solve(loads, fixed)
print(f"by the direct stiffness method:     {-stiffness.displacements[5] * 1e3:.3f} mm")
assert np.isclose(delta, -stiffness.displacements[5])

fig, ax = plt.subplots(figsize=(7, 3.5))
ax.bar([f"{i}-{j}" for i, j in members], contributions * 1e3)
ax.set(xlabel="member", ylabel="share of deflection [mm]", title="Where the midspan deflection comes from")
fig.tight_layout()

# %%
# The same idea for a beam
# ------------------------
#
# Midspan deflection of a simply supported beam under a uniform load: the
# unit load at midspan gives the triangle :math:`m = x/2` up to midspan.

L, w, I = 6.0, 10e3, 5e-5
x = np.linspace(0, L, 601)
M = w * x * (L - x) / 2
m = np.minimum(x, L - x) / 2
beam = unit_load_beam_deflection(x, M, m, E, I)
print(f"beam midspan deflection {beam * 1e3:.3f} mm, 5wL^4/384EI = {5 * w * L**4 / (384 * E * I) * 1e3:.3f} mm")
assert np.isclose(beam, 5 * w * L**4 / (384 * E * I))

plt.show()
