r"""
Hardy Cross's moment distribution (1930)
========================================

Continuous beams and building frames of the 1920s led to large systems of
simultaneous equations that no engineer wanted to solve by hand. In 1930
Hardy Cross proposed relaxation instead. Lock every joint against rotation
and apply the fixed-end moments. Then release one joint at a time: its
unbalanced moment is shared among the members meeting there in proportion
to their stiffness :math:`4EI/L`, and half of each share is *carried over*
to the far end of the member. Repeat until the joints balance
(:func:`~engineeringkit.structural_engineering.moment_distribution`). The
method fitted on a sheet of paper and was used for forty years. The example
balances a three-span beam built in at its left end, checks the result
against the plane-frame stiffness solution and plots how fast the
unbalance dies away.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import PlanarFrame, moment_distribution

# %%
# A three-span beam
# -----------------
#
# Joint 0 is built in; joints 1, 2 and 3 are on rollers. The middle span
# is stiffer and the last span unloaded.

spans = np.array([6.0, 8.0, 5.0])
EI = np.array([1.0, 2.0, 1.0]) * 40e6
w = np.array([20e3, 30e3, 0.0])
members = [(k, k + 1, spans[k], EI[k], w[k]) for k in range(3)]
result = moment_distribution(members, fixed_joints=(0,))
print(f"balanced in {result.iterations} cycles")
for k, (Mi, Mj) in enumerate(result.end_moments):
    print(f"span {k + 1}: end moments {Mi / 1e3:8.2f} and {Mj / 1e3:8.2f} kN m (clockwise positive)")
k0, k1 = 4 * EI[0] / spans[0], 4 * EI[1] / spans[1]
print(f"first cycle by hand: distribution factors at joint 1 are {k0 / (k0 + k1):.3f} and {k1 / (k0 + k1):.3f}")

# %%
# Check against the stiffness method
# ----------------------------------
#
# The plane-frame solver gives counterclockwise end moments, the negative
# of Cross's clockwise ones.

x = np.concatenate([[0], np.cumsum(spans)])
nodes = np.column_stack([x, np.zeros(4)])
frame = PlanarFrame(nodes, ((0, 1), (1, 2), (2, 3)), 1.0, 1e6, EI)
fixed_dofs = [0, 1, 2, 4, 7, 10]
stiffness = frame.solve(np.zeros(12), fixed_dofs, member_loads=w)
direct = -stiffness.member_end_forces[:, [2, 5]]
print(f"largest difference from the stiffness method: {np.max(abs(direct - result.end_moments)):.2e} N m")
assert np.allclose(result.end_moments, direct, rtol=1e-8, atol=1e-4)
for joint in (1, 2):
    assert np.isclose(result.end_moments[joint - 1, 1] + result.end_moments[joint, 0], 0, atol=1e-4)  # balanced

fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
result.plot(ax=axes[0])
axes[0].set_title("Unbalance at the joints")
for k in range(3):
    xs = np.linspace(0, spans[k], 50)
    Mi, Mj = result.end_moments[k, 0], -result.end_moments[k, 1]  # sagging positive
    sag = Mi + (Mj - Mi) * xs / spans[k] + w[k] * xs * (spans[k] - xs) / 2
    axes[1].plot(x[k] + xs, sag / 1e3, color="C0")
axes[1].axhline(0, color="gray", lw=0.5)
axes[1].invert_yaxis()
axes[1].set(xlabel="x [m]", ylabel="bending moment [kN m] (sagging down)", title="Bending moment diagram")
fig.tight_layout()

plt.show()
