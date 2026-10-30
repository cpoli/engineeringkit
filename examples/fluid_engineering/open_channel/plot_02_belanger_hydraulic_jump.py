r"""
Bélanger's hydraulic jump (1828)
================================

Water shooting fast and shallow from under a sluice gate cannot slow down
smoothly: it rises abruptly, in a turbulent roller, to a deeper, slower
stream. Jean-Baptiste Bélanger studied this hydraulic jump in his 1828
essay on steady flow in open channels. The rule used today, that the
pressure force plus the momentum flux is the same on both sides, gives the
sequent depth
:math:`y_2/y_1 = (\sqrt{1 + 8\mathrm{Fr}_1^2} - 1)/2`
(:func:`~engineeringkit.fluid_engineering.sequent_depth`). The energy
balance cannot be used across the jump, because the roller dissipates
energy (:func:`~engineeringkit.fluid_engineering.hydraulic_jump_head_loss`).
This example finds the jump below a sluice gate and how much of the
stream's energy it destroys.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import G, froude_number, hydraulic_jump_head_loss, sequent_depth, specific_energy

# %%
# The jump below a sluice gate
# ----------------------------
#
# A gate releases 4 m²/s per metre of width at a depth of 0.4 m.

q, y1 = 4.0, 0.4
v1 = q / y1
y2 = sequent_depth(y1, v1)
loss = hydraulic_jump_head_loss(y1, y2)
E1 = specific_energy(y1, q)
print(f"Fr1 = {froude_number(v1, y1):.2f}: depth rises from {y1} m to {y2:.2f} m (Fr2 = {froude_number(q / y2, y2):.2f})")
print(f"head lost in the jump {loss:.2f} m, {100 * loss / E1:.0f}% of the upstream specific energy")

# %%
# Momentum, not energy
# --------------------
#
# The momentum function :math:`M = y^2/2 + q^2/(gy)` is equal at the two
# depths, while the specific energy drops.

M1, M2 = y1**2 / 2 + q**2 / (G * y1), y2**2 / 2 + q**2 / (G * y2)
print(f"M1 = {M1:.3f} m², M2 = {M2:.3f} m²; E1 = {E1:.2f} m, E2 = {specific_energy(y2, q):.2f} m")

# %%
# Stronger jumps
# --------------
#
# Stilling basins below spillways are designed around the upstream Froude
# number: a weak undular jump near 1 loses little energy, a steady jump at
# 4.5 to 9 loses most of it.

Fr1 = np.linspace(1, 10, 100)
ratio = sequent_depth(1.0, Fr1 * np.sqrt(G))
loss_fraction = hydraulic_jump_head_loss(1.0, ratio) / specific_energy(1.0, Fr1 * np.sqrt(G))
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(Fr1, ratio, label="depth ratio $y_2/y_1$")
ax.plot(Fr1, 10 * loss_fraction, label="energy lost × 10")
ax.axvspan(4.5, 9, color="0.9", zorder=0, label="steady jump")
ax.set(xlabel="upstream Froude number $Fr_1$", title="Bélanger's hydraulic jump")
ax.legend()
fig.tight_layout()

assert np.isclose(M1, M2)
assert np.isclose(loss, E1 - specific_energy(y2, q))
assert froude_number(q / y2, y2) < 1 < froude_number(v1, y1)
assert np.isclose(ratio[0], 1.0) and loss_fraction[0] == 0

plt.show()
