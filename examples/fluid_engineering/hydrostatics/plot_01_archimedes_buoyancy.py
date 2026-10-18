r"""
Archimedes's principle (c. 250 BC)
==================================

In *On Floating Bodies* (c. 250 BC) Archimedes proved that a body immersed
in a fluid is pushed up by the weight of the fluid it displaces,
:math:`F_B = \rho g V` (:func:`~engineeringkit.fluid_engineering.buoyancy_force`).
Galileo's reconstruction of the story of King Hiero's crown uses that
principle: hang the crown and an equal mass of pure gold from a balance and
lower both into water. A crown alloyed with silver is bulkier, is pushed up
harder, and rises. This example weighs the crown in water, reads its silver
content back from the apparent weight, and floats an iceberg.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import G, buoyancy_force

# %%
# Weighing a crown in water
# -------------------------
#
# A 1 kg crown is a mixture of gold (19 300 kg/m³) and silver
# (10 490 kg/m³). Its volume, and so its buoyancy, grows with the silver
# fraction by mass :math:`x`.

RHO_GOLD, RHO_SILVER, RHO_WATER, MASS = 19300.0, 10490.0, 998.0, 1.0


def crown_volume(x):
    return MASS * ((1 - x) / RHO_GOLD + x / RHO_SILVER)


x = np.linspace(0, 1, 101)
apparent = MASS * G - buoyancy_force(RHO_WATER, crown_volume(x))
gold_bar = MASS * G - buoyancy_force(RHO_WATER, crown_volume(0.0))
print(f"pure gold in water weighs {gold_bar:.4f} N; a 30% silver crown {apparent[30]:.4f} N")

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(100 * x, apparent, label="crown")
ax.axhline(gold_bar, color="k", ls="--", label="equal mass of gold")
ax.set(xlabel="silver content [% by mass]", ylabel="apparent weight in water [N]", title="Hiero's crown")
ax.legend()
fig.tight_layout()

# %%
# Reading the alloy back
# ----------------------
#
# The apparent weight is linear in :math:`x`, so one weighing in air and one
# in water fix the silver fraction.

measured = apparent[30]
volume = (MASS * G - measured) / (RHO_WATER * G)
x_back = (volume / MASS - 1 / RHO_GOLD) / (1 / RHO_SILVER - 1 / RHO_GOLD)
print(f"silver fraction from the weighing: {x_back:.3f}")

# %%
# A floating iceberg
# ------------------
#
# A floating body displaces its own weight, so the submerged fraction is the
# density ratio: about 89% for ice (917 kg/m³) in sea water (1025 kg/m³).

rho_ice, rho_sea, V_ice = 917.0, 1025.0, 1.0e6
submerged = rho_ice * V_ice / rho_sea
print(f"submerged fraction of an iceberg: {submerged / V_ice:.3f}")

assert np.all(np.diff(apparent) < 0) and np.isclose(x_back, 0.30)
assert np.isclose(buoyancy_force(rho_sea, submerged), rho_ice * V_ice * G)

plt.show()
