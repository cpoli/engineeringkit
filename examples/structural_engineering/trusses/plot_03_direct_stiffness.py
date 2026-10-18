r"""
Turner, Clough, Martin and Topp: the direct stiffness method
============================================================

At Boeing in the early 1950s, M. J. Turner, R. W. Clough, H. C. Martin and
L. J. Topp needed the stiffness of swept jet wings, which hand methods
could not give. Their 1956 paper wrote down a stiffness matrix for each
element, added the matrices into one global matrix :math:`K` according to
which joints each element connects, and solved :math:`K_{ff}u_f = F_f` on a
computer. This is the *direct stiffness method*, the seed of the finite
element method (:class:`~engineeringkit.structural_engineering.PlanarTruss`).
It handles statically indeterminate trusses as easily as determinate ones,
but it fails without warning near a mechanism. The solver therefore
reports :math:`\operatorname{cond}(K_{ff})` and warns when it is large.
The example checks a redundant bridge truss against Clapeyron's energy
theorem, then flattens a joint until the warning fires.
"""

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import PlanarTruss, strain_energy_axial
from engineeringkit.structural_engineering.visualizers import plot_truss

# %%
# A redundant Warren truss
# ------------------------
#
# Four panels with both diagonals in the middle panels: Maxwell's count is
# +2, so statics alone cannot give the member forces.

nodes = np.array([[0, 0], [3, 0], [6, 0], [9, 0], [12, 0], [3, 3], [6, 3], [9, 3]], float)
members = (
    (0, 1), (1, 2), (2, 3), (3, 4),  # bottom chord
    (5, 6), (6, 7),  # top chord
    (0, 5), (4, 7), (1, 5), (2, 6), (3, 7),  # end posts and verticals
    (5, 2), (1, 6), (6, 3), (2, 7),  # crossed diagonals
)  # fmt: skip
E, A = 200e9, 2e-3
loads = np.zeros(16)
loads[[3, 5, 7]] = -80e3
fixed = [0, 1, 9]
result = PlanarTruss(nodes, members, E, A).solve(loads, fixed)
summary = result.summary()
print(f"largest displacement {1e3 * summary['max_abs_displacement']:.2f} mm, cond(K_ff) = {summary['condition_number']:.1f}")

# %%
# Checks: equilibrium and Clapeyron's theorem
# -------------------------------------------
#
# The reactions balance the loads, and the work of the loads,
# :math:`\tfrac12 F\cdot u`, equals the strain energy stored in the bars.

lengths = np.hypot(*(nodes[[j for _, j in members]] - nodes[[i for i, _ in members]]).T)
U = np.sum(strain_energy_axial(result.member_forces, lengths, A, E))
work = 0.5 * loads @ result.displacements
print(f"external work {work:.2f} J, strain energy {U:.2f} J")
assert np.isclose(work, U)
assert np.isclose(result.reactions.sum(), -loads.sum())

# %%
# Conditioning near a mechanism
# -----------------------------
#
# Two bars meeting at a joint carry a vertical load by their slope. As the
# joint is pushed down into line with the supports, the truss approaches a
# mechanism and :math:`\operatorname{cond}(K_{ff})` grows like
# :math:`1/\text{rise}^2`.

rises = np.logspace(-1, -7, 13)
conds, warned = [], []
for rise in rises:
    pair = PlanarTruss([[0, 0], [1, rise], [2, 0]], ((0, 1), (1, 2)), E, A)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        conds.append(pair.solve([0, 0, 0, -1e3, 0, 0], [0, 1, 4, 5]).condition_number)
    warned.append(any(issubclass(w.category, RuntimeWarning) for w in caught))
conds, warned = np.array(conds), np.array(warned)
print(f"first warning at a rise of {rises[warned][0]:.0e} m, cond(K_ff) = {conds[warned][0]:.1e}")
assert np.all(warned == (conds > 1e12))
assert np.allclose(conds * rises**2, conds[0] * rises[0] ** 2, rtol=0.05)

fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))

plot_truss(nodes, members, displacements=result.displacements, forces=result.member_forces, scale=200, ax=axes[0])
axes[0].set_title("Deflected shape (x200), tension red")
axes[1].loglog(rises, conds, "o-")
axes[1].loglog(rises[warned], conds[warned], "o", color="C3", label="RuntimeWarning")
axes[1].axhline(1e12, color="gray", ls="--")
axes[1].invert_xaxis()
axes[1].set(xlabel="rise of the loaded joint [m]", ylabel=r"cond($K_{ff}$)", title="Approaching a mechanism")
axes[1].legend()
fig.tight_layout()

plt.show()
