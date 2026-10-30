r"""
Maney's slope-deflection method (1915)
======================================

George A. Maney (1915), studying secondary stresses in rigid frames at the
University of Minnesota, wrote each member-end moment
in terms of the rotations of its two ends and the rotation of its chord:

.. math::

    M_{ij} = \frac{2EI}{L}\left(2\theta_i + \theta_j - 3\psi\right) + M^F_{ij}.

One equilibrium equation per free joint then gives the rotations
(:func:`~engineeringkit.structural_engineering.slope_deflection`). Hardy
Cross's moment distribution (1930) solves the same equations by
relaxation, and the stiffness method is their matrix form. Johnson's
*Modern Framed Structures* credits Otto Mohr with the idea before Maney
found it independently. The example analyses a three-span bridge girder,
compares it with moment distribution, and lets the middle pier settle.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import moment_distribution, slope_deflection

# %%
# Three spans, built in at the left abutment
# ------------------------------------------

EI, w = 4e8, 25e3
spans = [10.0, 14.0, 10.0]
members = [(k, k + 1, L, EI, w) for k, L in enumerate(spans)]
result = slope_deflection(members, fixed_joints=(0,))
cross = moment_distribution(members, fixed_joints=(0,))
print("joint rotations [mrad]:", np.round(result.rotations * 1e3, 3))
print("end moments [kN m]:\n", np.round(result.end_moments / 1e3, 1))
print(f"moment distribution needed {cross.iterations} cycles for the same answer")
assert np.allclose(result.end_moments, cross.end_moments, atol=1e-3)

# %%
# The second pier settles 20 mm
# -----------------------------
#
# The two spans meeting at joint 2 turn: the left chord clockwise
# (:math:`\psi = \Delta/L`), the right one counterclockwise.

delta = 0.02
psi = [0.0, delta / spans[1], -delta / spans[2]]
settled = slope_deflection(members, fixed_joints=(0,), chord_rotations=psi)
change = (settled.end_moments - result.end_moments) / 1e3
print("change in end moments from the settlement [kN m]:\n", np.round(change, 1))
for joint in (1, 2, 3):
    at = np.array([[i, j] for i, j, *_ in members]) == joint
    assert np.isclose(settled.end_moments[at].sum(), 0, atol=1e-6)

fig, ax = plt.subplots(figsize=(7, 3.5))
x0 = np.concatenate([[0], np.cumsum(spans)])
for k, L in enumerate(spans):
    x = np.linspace(0, L, 101)
    for moments, style, label in ((result.end_moments, "-", "level piers"), (settled.end_moments, "--", "pier 2 settled")):
        left, right = moments[k, 0], -moments[k, 1]  # sagging-positive moments at the span ends
        sag = left + (right - left) * x / L + w * x * (L - x) / 2
        ax.plot(x0[k] + x, sag / 1e3, style, color="C0" if style == "-" else "C3", label=label if k == 0 else None)
ax.axhline(0, color="gray", lw=0.5)
ax.set(xlabel="x [m]", ylabel="bending moment [kN m]", title="Three-span girder by slope-deflection")
ax.legend()
fig.tight_layout()

plt.show()
