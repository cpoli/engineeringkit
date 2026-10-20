r"""
Denavit and Hartenberg's link parameters
========================================

Jacques Denavit and Richard Hartenberg (1955) proposed a notation for any
chain of links joined by lower pairs. Each link frame is placed so that
going to the next takes only four numbers: a rotation :math:`\theta` and a
shift :math:`d` along the joint axis :math:`z`, then a shift :math:`a` and a
twist :math:`\alpha` along the common normal :math:`x`. The step is a single
4 by 4 homogeneous matrix
(:func:`~engineeringkit.mechanical_engineering.dh_transform`), and the pose
of a robot's hand is the product of one such matrix per joint
(:func:`~engineeringkit.mechanical_engineering.forward_kinematics`). The
example checks the product against hand-derived trigonometry for a planar
arm and a spatial articulated arm.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import dh_transform, forward_kinematics
from engineeringkit.mechanical_engineering.visualizers import plot_linkage

# %%
# A planar three-link arm
# -----------------------
#
# All joint axes are parallel (:math:`\alpha = 0`, :math:`d = 0`), so the
# hand position is the sum of link vectors at cumulative angles.

lengths = np.array([0.4, 0.3, 0.15])
q = np.radians([30.0, 45.0, -60.0])
frames = forward_kinematics(np.column_stack((lengths, np.zeros(3), np.zeros(3), q)))
expected = np.array([np.sum(lengths * np.cos(np.cumsum(q))), np.sum(lengths * np.sin(np.cumsum(q)))])
print("planar hand position:", frames[-1][:2, 3].round(6).tolist(), "expected", expected.round(6).tolist())
assert np.allclose(frames[-1][:2, 3], expected)
assert np.isclose(np.arctan2(frames[-1][1, 0], frames[-1][0, 0]), np.sum(q))

# %%
# A spatial articulated arm
# -------------------------
#
# A waist joint about the vertical (twist 90 degrees to the shoulder axis)
# carries a two-link planar arm, as in the PUMA's first three joints. The
# hand lies at radius :math:`\rho = a_2\cos q_2 + a_3\cos(q_2 + q_3)` from
# the waist axis, height :math:`d_1 + a_2\sin q_2 + a_3\sin(q_2 + q_3)`.

d1, a2, a3 = 0.5, 0.45, 0.4
for q1, q2, q3 in np.radians([[0, 30, -45], [60, 75, -120], [-135, 10, 30]]):
    table = [[0, np.pi / 2, d1, q1], [a2, 0, 0, q2], [a3, 0, 0, q3]]
    hand = forward_kinematics(table)[-1]
    rho = a2 * np.cos(q2) + a3 * np.cos(q2 + q3)
    expected = [rho * np.cos(q1), rho * np.sin(q1), d1 + a2 * np.sin(q2) + a3 * np.sin(q2 + q3)]
    print(f"q = {np.degrees([q1, q2, q3]).round(0).tolist()}: hand {hand[:3, 3].round(4).tolist()}")
    assert np.allclose(hand[:3, 3], expected)
    assert np.allclose(hand[:3, :3] @ hand[:3, :3].T, np.eye(3)) and np.isclose(np.linalg.det(hand[:3, :3]), 1)

# %%
# The inverse step
# ----------------
#
# Each transform is a rigid motion, so its inverse is
# :math:`[R^T,\ -R^T p]` and the product of the two is the identity.

T = dh_transform(0.3, np.pi / 3, 0.1, 0.7)
inverse = np.eye(4)
inverse[:3, :3], inverse[:3, 3] = T[:3, :3].T, -T[:3, :3].T @ T[:3, 3]
assert np.allclose(T @ inverse, np.eye(4))

fig, ax = plt.subplots(figsize=(5.5, 4.5))
for angle in np.radians([0, 40, 80, 120]):
    pose = q + np.array([angle, -angle / 2, angle / 3])
    joints = forward_kinematics(np.column_stack((lengths, np.zeros(3), np.zeros(3), pose)))[:, :2, 3]
    plot_linkage(joints, ax=ax, alpha=0.8, label=f"first joint +{np.degrees(angle):.0f} deg")
ax.legend(fontsize=8)
ax.set(title="Planar arm poses from DH frames")
fig.tight_layout()

plt.show()
