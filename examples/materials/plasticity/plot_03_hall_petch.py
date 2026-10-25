r"""
Hall and Petch: grain-size strengthening (1951-1953)
====================================================

Eric Hall (1951) measured the lower yield stress of mild steels with
different grain sizes, and Norman Petch (1953) the cleavage strength of
iron. Both found the strength rising linearly with :math:`d^{-1/2}`,
:math:`\sigma_y = \sigma_0 + k\,d^{-1/2}`
(:func:`~engineeringkit.materials.hall_petch`). Dislocations pile up against
grain boundaries, and the stress concentrated at the head of a pile-up of
length :math:`d` grows as :math:`\sqrt d`; slip crosses into the next grain
when it reaches a fixed value. Grain refinement is the only strengthening
mechanism that also makes steel tougher.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import hall_petch

# %%
# Fitting the Hall-Petch line
# ---------------------------
#
# Yield strengths of a mild steel, :math:`\sigma_0 = 70` MPa and
# :math:`k = 0.74` MPa :math:`\sqrt{m}`, measured with 5 MPa scatter at
# grain sizes from 2 to 200 micrometres. A straight-line fit in
# :math:`d^{-1/2}` recovers both constants.

sigma0, k = 70e6, 0.74e6
d = np.geomspace(2e-6, 200e-6, 12)
rng = np.random.default_rng(2)
measured = hall_petch(sigma0, k, d) + rng.normal(0, 5e6, d.size)
k_fit, sigma0_fit = np.polyfit(d**-0.5, measured, 1)
print(f"fitted sigma_0 = {sigma0_fit / 1e6:.1f} MPa, k = {k_fit / 1e6:.3f} MPa sqrt(m)")
assert abs(k_fit / k - 1) < 0.05 and abs(sigma0_fit - sigma0) < 10e6

# %%
# How much refinement buys
# ------------------------
#
# Going from a 50 um to a 5 um grain more than doubles the yield strength
# of this steel. Below about 20 nm the trend reverses as grain boundaries
# start to slide, and the relation is not used there.

coarse, fine = hall_petch(sigma0, k, 50e-6), hall_petch(sigma0, k, 5e-6)
print(f"50 um: {coarse / 1e6:.0f} MPa; 5 um: {fine / 1e6:.0f} MPa")
assert fine > 2 * coarse

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
x = np.linspace(0, d.min() ** -0.5, 2)
ax1.plot(d**-0.5 / 1e3, measured / 1e6, "o", label="tests")
ax1.plot(x / 1e3, (sigma0_fit + k_fit * x) / 1e6, label="fit")
ax1.set(xlabel="d^-1/2 [mm^-1/2]", ylabel="yield strength [MPa]", title="Hall-Petch plot")
ax1.legend()
grains = np.geomspace(1e-6, 1e-3, 100)
ax2.semilogx(1e6 * grains, hall_petch(sigma0, k, grains) / 1e6)
ax2.set(xlabel="grain size [um]", ylabel="yield strength [MPa]", title="Mild steel")
fig.tight_layout()

plt.show()
