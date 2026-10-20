r"""
Whipple's method of joints
==========================

American railway bridges of the 1840s were built from rules of thumb, and
some fell down. Squire Whipple, in *A Work on Bridge Building* (1847),
showed that the force in every member of a pin-jointed truss follows from
statics alone. Each joint is in equilibrium, so its two equations,
:math:`\sum F_x = 0` and :math:`\sum F_y = 0`, give the forces in two
members. Start at a support and work across the truss, two unknowns at a
time (:func:`~engineeringkit.structural_engineering.method_of_joints`). The
example solves a six-panel Pratt truss under a train of panel loads and
checks the chord forces against the method of sections: at any vertical
cut, the chord force is the bending moment divided by the truss depth.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import influence_line, method_of_joints
from engineeringkit.structural_engineering.visualizers import plot_truss

# %%
# A Pratt truss with six panels
# -----------------------------
#
# Joints 0-6 run along the bottom chord, joints 7-11 along the top. The
# left end is pinned and the right end is on a roller.

panels, a, depth, P = 6, 4.0, 4.0, 50e3
bottom = [(a * k, 0.0) for k in range(panels + 1)]
top = [(a * k, depth) for k in range(1, panels)]
nodes = np.array(bottom + top)
T = {k: panels + k for k in range(1, panels)}  # top joint above bottom joint k
members = [(k, k + 1) for k in range(panels)]  # bottom chord
members += [(T[k], T[k + 1]) for k in range(1, panels - 1)]  # top chord
members += [(0, T[1]), (panels, T[panels - 1])]  # end posts
members += [(k, T[k]) for k in range(1, panels)]  # verticals
members += [(T[k], k + 1) for k in range(1, panels // 2)] + [(T[k], k - 1) for k in range(panels // 2 + 1, panels)]  # diagonals

loads = np.zeros(2 * len(nodes))
loads[2 * np.arange(1, panels) + 1] = -P  # a load at every interior bottom joint
fixed = [0, 1, 2 * panels + 1]
result = method_of_joints(nodes, members, loads, fixed)
left, right = result.reactions[1], result.reactions[2 * panels + 1]
print(f"{len(members)} members, {len(nodes)} joints; reactions {left / 1e3:.0f} and {right / 1e3:.0f} kN")
print("joints used in the order", result.order)
summary = result.summary()
print(f"{summary['tension_members']} members in tension, {summary['compression_members']} in compression")

# %%
# Check: the method of sections
# -----------------------------
#
# Cut the truss through the panel next to top joint :math:`k`, on the side
# its diagonal runs to. The top chord and the diagonal both pass through
# that joint, so moments about it give the bottom-chord tension as the
# bending moment there divided by the depth. At midspan, moments about the
# bottom joint give the top-chord compression the same way. The moment is
# the sum of the loads times the moment influence line at the section.

L = panels * a
for k in range(1, panels):
    moment = sum(P * influence_line(L, a * j, "moment", at=a * k) for j in range(1, panels))
    if k == panels // 2:
        name, chord = "top chord", -result.member_forces[panels + k - 2]
    else:
        name, chord = "bottom chord", result.member_forces[k if k < panels // 2 else k - 1]
    print(f"x = {a * k:4.0f} m: M/h = {moment / depth / 1e3:6.1f} kN, {name} {abs(chord) / 1e3:6.1f} kN")
    assert np.isclose(chord, moment / depth)
assert np.isclose(left + right, (panels - 1) * P)

fig, ax = plt.subplots(figsize=(8, 3))
plot_truss(nodes, members, forces=result.member_forces, ax=ax)
for (i, j), force in zip(members, result.member_forces):
    mx, my = nodes[[i, j]].mean(axis=0)
    ax.annotate(f"{force / 1e3:.0f}", (mx, my), fontsize=7, ha="center", va="center", backgroundcolor="w")
ax.set_title("Pratt truss: tension red, compression blue [kN]")
fig.tight_layout()

plt.show()
