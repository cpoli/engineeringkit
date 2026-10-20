r"""
Colebrook's turbulent friction factor
=====================================

Commercial pipes do not behave like Nikuradse's sand: their roughness is
uneven, and friction leaves the smooth law gradually rather than with
Nikuradse's dip. Cyril Colebrook, after experiments with Cedric White on
pipes coated with mixed grains, added the smooth and rough laws inside one
logarithm in 1939,

.. math::

    \frac{1}{\sqrt f} = -2\log_{10}\left(\frac{\varepsilon/D}{3.7} + \frac{2.51}{Re\sqrt f}\right),

(:func:`~engineeringkit.fluid_engineering.colebrook_friction_factor`). The
equation is implicit in :math:`f`. This example solves it by the
fixed-point iteration engineers used by hand and compares it with the two
limits it joins.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import (
    colebrook_friction_factor,
    nikuradse_friction_factor,
    smooth_pipe_friction_factor,
)

# %%
# Fixed-point iteration
# ---------------------
#
# With :math:`x = 1/\sqrt f` the equation reads
# :math:`x = -2\log_{10}(\varepsilon/3.7D + 2.51x/Re)`. The right side
# changes slowly with :math:`x`, so substitution converges in a few steps.

Re, rel = 2e5, 1e-3
x = [3.0]
for _ in range(8):
    x.append(-2 * np.log10(rel / 3.7 + 2.51 * x[-1] / Re))
f_iter = 1 / np.array(x) ** 2
f_ref = colebrook_friction_factor(Re, rel)
for k, fk in enumerate(f_iter[:5]):
    print(f"step {k}: f = {fk:.6f}")
print(f"bracketed root: f = {f_ref:.6f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.8))
ax1.semilogy(np.abs(f_iter - f_ref) + 1e-17, marker="o")
ax1.set(xlabel="iteration", ylabel="|f - f*|", title="Substitution converges fast")

# %%
# A smooth bridge between two laws
# --------------------------------
#
# Colebrook's curve sits above both the smooth and the fully rough laws and
# approaches each at its end of the range. The largest excess is in the
# transition, where both mechanisms count.

Re_range = np.geomspace(4e3, 1e8, 200)
cb = colebrook_friction_factor(Re_range, rel)
smooth = smooth_pipe_friction_factor(Re_range)
rough = nikuradse_friction_factor(rel)
ax2.loglog(Re_range, cb, label="Colebrook, ε/D = 0.001")
ax2.loglog(Re_range, smooth, "k--", label="smooth")
ax2.axhline(rough, color="k", ls=":", label="fully rough")
ax2.set(xlabel="Reynolds number", ylabel="friction factor", title="Colebrook joins the limits")
ax2.legend(fontsize="small")
fig.tight_layout()

excess = cb / np.maximum(smooth, rough) - 1
print(f"largest excess over the larger limit: {excess.max():.0%} at Re = {Re_range[excess.argmax()]:.1e}")

assert abs(f_iter[-1] - f_ref) < 1e-10 * f_ref
assert np.all(cb >= np.maximum(smooth, rough) * (1 - 2e-3))
assert abs(cb[0] / smooth[0] - 1) < 0.12 and abs(cb[-1] / rough - 1) < 0.01

plt.show()
