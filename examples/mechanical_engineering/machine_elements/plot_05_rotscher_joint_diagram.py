r"""
Rötscher's bolted-joint diagram (1927)
======================================

Nineteenth-century designers sized a bolt by dividing the load by its
root area. Friedrich Rötscher's *Die Maschinenelemente* (1927) drew the
joint diagram that showed why a preloaded bolt does far better: bolt and
clamped members are two springs, one stretched and one squeezed by the
preload :math:`F_i`, and an external load :math:`P` stretches the bolt
further only by as much as it relieves the members. The bolt sees only the
fraction :math:`C = k_b/(k_b + k_m)` of the load
(:func:`~engineeringkit.mechanical_engineering.bolted_joint`) until the
members separate. The members' stiffness comes from cones of compressed
material spreading from under the head and nut
(:func:`~engineeringkit.mechanical_engineering.member_stiffness`), the
bolt's from its shank and threads in series
(:func:`~engineeringkit.mechanical_engineering.bolt_stiffness`). Rötscher's
1927 book could not be consulted; the date and attribution follow J. H.
Bickford, *An Introduction to the Design and Behavior of Bolted Joints*,
and the Bolt Council's history of joint design.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

from engineeringkit.mechanical_engineering import bolt_stiffness, bolted_joint, member_stiffness, preload_torque

# %%
# Stiffness of an M12 bolt and two steel plates
# ---------------------------------------------
#
# The frustum formula is the exact compliance of two 30-degree cones; check
# it by integrating :math:`dx/EA(x)`.

E, d, grip = 207e9, 0.012, 0.04
At, Ad = 84.3e-6, np.pi * d**2 / 4  # M12 x 1.75 tensile-stress area
lt, ld = 0.018, grip - 0.018
kb = bolt_stiffness(E, Ad, ld, At, lt)
km = member_stiffness(E, d, grip)
cone = 2 * quad(lambda x: 1 / (E * np.pi * ((0.75 * d + x * np.tan(np.pi / 6)) ** 2 - (d / 2) ** 2)), 0, grip / 2)[0]
print(f"kb = {kb / 1e9:.3f} GN/m, km = {km / 1e9:.3f} GN/m (cones integrated: {1 / cone / 1e9:.3f}), C = {kb / (kb + km):.3f}")
assert np.isclose(km, 1 / cone)

# %%
# Preload and the load factor
# ---------------------------
#
# Property class 8.8: proof strength 600 MPa. Preload to 75% of the proof
# load for a reusable joint, then apply a 15 kN separating load.

proof = 600e6 * At
Fi, P = 0.75 * proof, 15e3
joint = bolted_joint(kb, km, Fi, P)
print({key: round(value, 3) for key, value in joint.summary().items()})
print(f"load factor against proof: {joint.load_factor(proof):.2f}; wrench torque {preload_torque(Fi, d):.0f} N m (K = 0.2)")
# The bolt load grows by only C P, so the load factor is much larger than (proof - Fi)/P.
assert np.isclose(joint.bolt_load - Fi, joint.joint_constant * P)
assert joint.load_factor(proof) > (proof - Fi) / P
assert joint.separation_load > P

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
joint.plot(ax=ax1)
ax1.set(title="Joint diagram")
loads = np.linspace(1.0, 2 * joint.separation_load, 300)
ax2.plot(loads / 1e3, [bolted_joint(kb, km, Fi, p).bolt_load / 1e3 for p in loads], label="bolt tension")
ax2.plot(loads / 1e3, [bolted_joint(kb, km, Fi, p).member_load / 1e3 for p in loads], label="member clamp")
ax2.axvline(joint.separation_load / 1e3, color="gray", ls="--", label="separation")
ax2.set(xlabel="external load [kN]", ylabel="force [kN]")
ax2.legend()
fig.tight_layout()

plt.show()
