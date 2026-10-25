r"""
Maxwell's rigidity count and reciprocal theorem (1864)
======================================================

James Clerk Maxwell's short 1864 paper "On the calculation of the
equilibrium and stiffness of frames" made two contributions. First, a
counting rule: a plane frame of :math:`j` joints needs :math:`2j - 3`
bars to be stiff, or, counting the :math:`r` support reactions,
:math:`m + r = 2j` for a statically determinate truss
(:func:`~engineeringkit.structural_engineering.truss_determinacy`). The
count is necessary but not sufficient, since bars can be placed where they
are not needed. Second, the reciprocal theorem: the deflection at joint
:math:`i` caused by a unit load at joint :math:`j` equals the deflection at
:math:`j` caused by a unit load at :math:`i`
(:func:`~engineeringkit.structural_engineering.maxwell_betti_reciprocity`).
Betti generalised it in 1872. The example counts four trusses, shows one
that passes the count and is still a mechanism, and checks reciprocity on
a statically indeterminate truss.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import PlanarTruss, maxwell_betti_reciprocity, truss_determinacy
from engineeringkit.structural_engineering.visualizers import plot_truss

# %%
# Counting bars
# -------------

square = np.array([[0, 0], [2, 0], [2, 2], [0, 2]], float)
sides = [(0, 1), (1, 2), (2, 3), (3, 0)]
cases = {
    "unbraced square": (sides, [0, 1, 3]),
    "braced square": (sides + [(0, 2)], [0, 1, 3]),
    "double-braced square": (sides + [(0, 2), (1, 3)], [0, 1, 3]),
    "unbraced square on two pins": (sides, [0, 1, 2, 3]),
}
fig, axes = plt.subplots(1, 4, figsize=(12, 3))
for ax, (name, (members, fixed)) in zip(axes, cases.items()):
    count = truss_determinacy(len(square), len(members), len(fixed))
    try:
        PlanarTruss(square, tuple(members), 200e9, 1e-3).solve(np.zeros(8), fixed)
        verdict = "stiff"
    except ValueError:
        verdict = "mechanism"
    print(f"{name:28s} m + r - 2j = {count:+d}: {verdict}")
    plot_truss(square, members, ax=ax)
    ax.plot(*square[sorted({d // 2 for d in fixed})].T, "^", color="C2", ms=10)
    ax.set_title(f"{name}\ncount {count:+d}, {verdict}", fontsize=9)
    ax.set_axis_off()
    if name == "unbraced square on two pins":
        assert count == 0 and verdict == "mechanism"  # the count is not sufficient
fig.tight_layout()

# %%
# The reciprocal theorem
# ----------------------
#
# The double-braced square on a pin and a roller has one redundant bar.
# Each column of the flexibility matrix is the displacement field under a
# unit load at one degree of freedom; Maxwell's theorem says the matrix is
# symmetric.

members, fixed = cases["double-braced square"]
truss = PlanarTruss(square, tuple(members), 200e9, 1e-4)
free = [dof for dof in range(8) if dof not in fixed]
flex = np.empty((len(free), len(free)))
for col, dof in enumerate(free):
    unit = np.zeros(8)
    unit[dof] = 1.0
    flex[:, col] = truss.solve(unit, fixed).displacements[free]
i, j = free.index(4), free.index(7)  # joint 2 horizontally, joint 3 vertically
print(f"delta_ij = {flex[i, j]:.4e} m/N, delta_ji = {flex[j, i]:.4e} m/N")
assert maxwell_betti_reciprocity(flex[i, j], flex[j, i])
assert np.allclose(flex, flex.T, rtol=1e-12, atol=1e-20)
print(f"largest asymmetry of the flexibility matrix: {np.max(abs(flex - flex.T)):.1e} m/N")

plt.show()
